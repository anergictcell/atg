use std::cmp::Ordering;
use std::cmp::{max, min};
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenomicRelation {
    Match,
    Upstream,
    Downstream,
    Overlaps,
    Inside,
    Left,
    Right,
}

impl fmt::Display for GenomicRelation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                GenomicRelation::Match => "match",
                GenomicRelation::Upstream => "upstream",
                GenomicRelation::Downstream => "downstream",
                GenomicRelation::Overlaps => "overlaps",
                GenomicRelation::Inside => "inside",
                GenomicRelation::Left => "left",
                GenomicRelation::Right => "light",
            }
        )
    }
}

/// Returns the intersection coordinates between two genomic features
///
/// ```text
/// a:       --AAA-- --AAA-- --AAA-- --AAA-- --AAA-- --AAA-- --AAA--
/// b:       BB----- -----BB -BBBBB- --BBB-- ---B--- -BB---- ----BB-
/// Returns: ------- ------- --XXX-- --XXX-- ---X--- --X---- ----X--
/// ```
///
/// ```rust
/// use atg::utils::intersect;
/// let a = (&3, &5);
/// let b = (&2, &6);
/// let c = intersect(a, b).unwrap();
/// assert_eq!(c, (3, 5));
/// ```
pub fn intersect(a: (&u32, &u32), b: (&u32, &u32)) -> Option<(u32, u32)> {
    if a.0 <= b.1 && a.1 >= b.0 {
        // Some overlap between both features
        Some((max(*a.0, *b.0), min(*a.1, *b.1)))
    } else {
        None
    }
}

/// Returns the union coordinates between two overlapping genomic features
///
/// ```text
/// a:       --AAA-- --AAA-- --AAA-- --AAA-- --AAA-- --AAA-- --AAA--
/// b:       BB----- -----BB -BBBBB- --BBB-- ---B--- -BB---- ----BB-
/// Returns: ------- ------- -XXXXX- --XXX-- --XXX-- -XXXX-- --XXX--
/// ```
///
/// ```rust
/// use atg::utils::union;
/// let a = (&3, &5);
/// let b = (&2, &3);
/// let c = union(a, b).unwrap();
/// assert_eq!(c, (2, 5));
/// ```
pub fn union(a: (&u32, &u32), b: (&u32, &u32)) -> Option<(u32, u32)> {
    if a.0 <= b.1 && a.1 >= b.0 {
        // Some overlap between both features
        Some((min(*a.0, *b.0), max(*a.1, *b.1)))
    } else {
        None
    }
}

/// Returns genomic features that intersect with a, but not with b
/// If feature b overlaps a, returns None
///
/// ```text
/// a:       --AAA-- --AAA-- --AAA-- --AAA-- --AAA-- --AAA-- --AAA--
/// b:       BB----- -----BB -BBBBB- --BBB-- ---B--- -BB---- ----BB-
/// Returns: --XXX-- --XXX-- ------- ------- --X-X-- ---XX-- --XX---
/// ```
///
/// ```rust
/// use atg::utils::subtract;
/// let a = (&3, &5);
/// let b = (&2, &3);
/// let c = subtract(a, b);
/// assert_eq!(c, vec![(4, 5)]);
/// ```
pub fn subtract(a: (&u32, &u32), b: (&u32, &u32)) -> Vec<(u32, u32)> {
    match relation(a, b) {
        GenomicRelation::Upstream => vec![(*a.0, *a.1)],
        GenomicRelation::Downstream => vec![(*a.0, *a.1)],
        GenomicRelation::Match => vec![],
        GenomicRelation::Overlaps => vec![(*a.0, *b.0 - 1), (*b.1 + 1, *a.1)],
        GenomicRelation::Inside => vec![],
        GenomicRelation::Left => vec![(*a.0, *b.0 - 1)],
        GenomicRelation::Right => vec![(*b.1 + 1, *a.1)],
    }
}

/// Returns the genomic relation of two features to each other
///
/// ```text
/// --AAA--
/// BB-----
/// Downstream
///
/// --AAA--
/// -----BB
/// Upstream
///
/// --AAA--
/// -BBBBB-
/// Inside
///
/// --AAA--
/// --BBB--
/// Match
///
/// --AAA--
/// ---B---
/// Overlaps
///
/// --AAA--
/// -BB----
/// Right
///
/// --AAA--
/// ----BB-
/// Left
/// ```
///
/// ```rust
/// use atg::utils::relation;
/// use atg::utils::GenomicRelation;
/// let a = (&3, &5);
/// let b = (&2, &3);
/// let c = relation(a, b);
/// assert_eq!(c, GenomicRelation::Right);
/// ```
pub fn relation(a: (&u32, &u32), b: (&u32, &u32)) -> GenomicRelation {
    if a.1 < b.0 {
        return GenomicRelation::Upstream;
    }
    if a.0 > b.1 {
        return GenomicRelation::Downstream;
    }
    match (a.0.cmp(b.0), a.1.cmp(b.1)) {
        (Ordering::Equal, Ordering::Equal) => GenomicRelation::Match,
        (start, end) if start.is_le() && end.is_ge() => GenomicRelation::Overlaps,
        (start, end) if start.is_ge() && end.is_le() => GenomicRelation::Inside,
        (start, end) if start.is_le() && end.is_le() => GenomicRelation::Left,
        (start, end) if start.is_ge() && end.is_ge() => GenomicRelation::Right,
        _ => panic!("Unable to determine relation of two genomic features"),
    }
}

/// Calculates the intersect between an exon and the coding sequence
///
/// # Examples
///
/// ```rust
/// use atg::utils::exon_cds_overlap;
/// let cds_start = 4;
/// let cds_end = 7;
/// let exon_start = 1;
/// let exon_end = 10;
/// assert_eq!(
///     exon_cds_overlap(&exon_start, &exon_end, &cds_start, &cds_end),
///     (Some(4), Some(7))
/// );
///
/// ```
pub fn exon_cds_overlap(
    exon_start: &u32,
    exon_end: &u32,
    cds_start: &u32,
    cds_end: &u32,
) -> (Option<u32>, Option<u32>) {
    match intersect((exon_start, exon_end), (cds_start, cds_end)) {
        Some((x, y)) => (Some(x), Some(y)),
        _ => (None, None),
    }
}

#[cfg(test)]
mod test_intersect {
    use super::intersect;

    #[test]
    fn test_no_overlap() {
        assert_eq!(intersect((&3, &5), (&1, &1)), None);
        assert_eq!(intersect((&3, &5), (&1, &2)), None);
        assert_eq!(intersect((&3, &5), (&2, &2)), None);
        assert_eq!(intersect((&3, &5), (&6, &6)), None);
        assert_eq!(intersect((&3, &5), (&6, &7)), None);
    }

    #[test]
    fn test_a_fully_in_b() {
        assert_eq!(intersect((&3, &5), (&1, &5)), Some((3, 5)));
        assert_eq!(intersect((&3, &5), (&1, &6)), Some((3, 5)));
        assert_eq!(intersect((&3, &5), (&3, &5)), Some((3, 5)));
        assert_eq!(intersect((&3, &5), (&3, &6)), Some((3, 5)));
    }

    #[test]
    fn test_b_fully_in_a() {
        assert_eq!(intersect((&1, &5), (&3, &5)), Some((3, 5)));
        assert_eq!(intersect((&1, &6), (&3, &5)), Some((3, 5)));
        assert_eq!(intersect((&3, &5), (&3, &5)), Some((3, 5)));
        assert_eq!(intersect((&3, &6), (&3, &5)), Some((3, 5)));
    }

    #[test]
    fn test_left_overlap() {
        assert_eq!(intersect((&3, &5), (&1, &3)), Some((3, 3)));
        assert_eq!(intersect((&3, &5), (&1, &4)), Some((3, 4)));
    }

    #[test]
    fn test_right_overlap() {
        assert_eq!(intersect((&3, &5), (&4, &6)), Some((4, 5)));
        assert_eq!(intersect((&3, &5), (&5, &6)), Some((5, 5)));
        assert_eq!(intersect((&3, &5), (&5, &7)), Some((5, 5)));
    }
}

#[cfg(test)]
mod test_relation {
    use super::relation;
    use super::GenomicRelation;

    #[test]
    fn test_upstream() {
        assert_eq!(relation((&3, &5), (&6, &7)), GenomicRelation::Upstream);
        assert_eq!(relation((&1, &2), (&3, &5)), GenomicRelation::Upstream);
        assert_eq!(relation((&1, &1), (&2, &4)), GenomicRelation::Upstream);
    }

    #[test]
    fn test_downstream() {
        assert_eq!(relation((&3, &5), (&1, &2)), GenomicRelation::Downstream);
        assert_eq!(relation((&6, &7), (&3, &5)), GenomicRelation::Downstream);
        assert_eq!(relation((&2, &4), (&1, &1)), GenomicRelation::Downstream);
    }

    #[test]
    fn test_inside() {
        assert_eq!(relation((&3, &5), (&2, &5)), GenomicRelation::Inside);
        assert_eq!(relation((&3, &5), (&2, &6)), GenomicRelation::Inside);
        assert_eq!(relation((&3, &5), (&3, &6)), GenomicRelation::Inside);
        assert_eq!(relation((&3, &3), (&3, &6)), GenomicRelation::Inside);
        assert_eq!(relation((&4, &4), (&3, &6)), GenomicRelation::Inside);
        assert_eq!(relation((&6, &6), (&3, &6)), GenomicRelation::Inside);
    }

    #[test]
    fn test_match() {
        assert_eq!(relation((&3, &5), (&3, &5)), GenomicRelation::Match);
        assert_eq!(relation((&3, &3), (&3, &3)), GenomicRelation::Match);
    }

    #[test]
    fn test_overlaps() {
        assert_eq!(relation((&2, &5), (&3, &5)), GenomicRelation::Overlaps);
        assert_eq!(relation((&2, &6), (&3, &5)), GenomicRelation::Overlaps);
        assert_eq!(relation((&3, &6), (&3, &5)), GenomicRelation::Overlaps);
        assert_eq!(relation((&3, &6), (&3, &3)), GenomicRelation::Overlaps);
        assert_eq!(relation((&3, &6), (&4, &4)), GenomicRelation::Overlaps);
        assert_eq!(relation((&3, &6), (&6, &6)), GenomicRelation::Overlaps);
    }

    #[test]
    fn test_right() {
        assert_eq!(relation((&3, &5), (&2, &3)), GenomicRelation::Right);
        assert_eq!(relation((&3, &5), (&2, &4)), GenomicRelation::Right);
        assert_eq!(relation((&3, &4), (&2, &3)), GenomicRelation::Right);
    }

    #[test]
    fn test_left() {
        assert_eq!(relation((&3, &5), (&5, &6)), GenomicRelation::Left);
        assert_eq!(relation((&1, &2), (&2, &4)), GenomicRelation::Left);
        assert_eq!(relation((&2, &3), (&3, &5)), GenomicRelation::Left);
        assert_eq!(relation((&2, &4), (&3, &5)), GenomicRelation::Left);
        assert_eq!(relation((&2, &3), (&3, &4)), GenomicRelation::Left);
    }
}

#[cfg(test)]
mod test_subtraction {
    use super::subtract;

    #[test]
    fn test_case_1() {
        assert_eq!(subtract((&3, &5), (&1, &2)), vec![(3, 5)]);
    }

    #[test]
    fn test_case_2() {
        assert_eq!(subtract((&3, &5), (&6, &7)), vec![(3, 5)]);
    }

    #[test]
    fn test_case_3() {
        assert_eq!(subtract((&3, &5), (&2, &6)), vec![]);
    }

    #[test]
    fn test_case_4() {
        assert_eq!(subtract((&3, &5), (&3, &5)), vec![]);
    }

    #[test]
    fn test_case_5() {
        assert_eq!(subtract((&3, &5), (&4, &4)), vec![(3, 3), (5, 5)]);
    }

    #[test]
    fn test_case_6() {
        assert_eq!(subtract((&3, &5), (&2, &3)), vec![(4, 5)]);
    }

    #[test]
    fn test_case_7() {
        assert_eq!(subtract((&3, &5), (&5, &6)), vec![(3, 4)]);
    }
}

#[cfg(test)]
mod test_overlap {
    use super::*;

    #[test]
    fn test_cds_overlap() {
        assert_eq!(exon_cds_overlap(&5, &8, &1, &10), (Some(5), Some(8)));
        assert_eq!(exon_cds_overlap(&5, &8, &5, &10), (Some(5), Some(8)));
        assert_eq!(exon_cds_overlap(&5, &8, &7, &10), (Some(7), Some(8)));
        assert_eq!(exon_cds_overlap(&5, &8, &1, &7), (Some(5), Some(7)));
        assert_eq!(exon_cds_overlap(&5, &8, &1, &8), (Some(5), Some(8)));
        assert_eq!(exon_cds_overlap(&5, &8, &7, &8), (Some(7), Some(8)));
        assert_eq!(exon_cds_overlap(&5, &8, &6, &7), (Some(6), Some(7)));
        assert_eq!(exon_cds_overlap(&5, &8, &2, &4), (None, None));
        assert_eq!(exon_cds_overlap(&5, &8, &9, &11), (None, None));
        assert_eq!(exon_cds_overlap(&5, &8, &8, &8), (Some(8), Some(8)));
        assert_eq!(exon_cds_overlap(&5, &8, &5, &5), (Some(5), Some(5)));
        assert_eq!(exon_cds_overlap(&5, &8, &6, &6), (Some(6), Some(6)));
    }
}
