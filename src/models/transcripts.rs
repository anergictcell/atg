use std::collections::HashMap;

use crate::models::Transcript;

/// A convinience wrapper to handle  large amounts of `Transcript`s
///
/// It allows fast lookup operation by gene or transcript name.
///
/// # Examples
///
/// ```rust
/// use atg::models::{TranscriptBuilder, Transcripts};
///
/// let mut transcripts = Transcripts::new();
/// assert_eq!(transcripts.len(), 0);
///
/// transcripts.push(TranscriptBuilder::new()
///     .name("NM_001203247.2")
///     .chrom("chr7")
///     .gene("EZH2")
///     .strand(atg::models::Strand::Minus)
///     .build()
///     .unwrap()
/// );
/// assert_eq!(transcripts.len(), 1);
///
/// assert!(transcripts.by_name("NM_001203247.2").is_some());
/// assert!(transcripts.by_gene("EZH2").is_some());
/// assert_eq!(transcripts.by_gene("EZH2").unwrap().len(), 1);
///
/// assert!(transcripts.by_name("Foo").is_none());
/// assert!(transcripts.by_gene("Bar").is_none());
/// ```
pub struct Transcripts {
    list: Vec<Transcript>,
    name: HashMap<String, usize>,
    gene: HashMap<String, Vec<usize>>,
}

impl Transcripts {
    pub fn new() -> Self {
        Self {
            list: vec![],
            name: HashMap::new(),
            gene: HashMap::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            list: Vec::with_capacity(capacity),
            name: HashMap::with_capacity(capacity),
            gene: HashMap::with_capacity(capacity),
        }
    }

    /// Retrieve a [`Transcript`] by its name / transcript-id
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use atg::models::{TranscriptBuilder, Transcripts};
    /// # let mut transcripts = Transcripts::new();
    /// # transcripts.push(TranscriptBuilder::new()
    /// #     .name("NM_001203247.2")
    /// #     .chrom("chr7")
    /// #     .gene("EZH2")
    /// #     .strand(atg::models::Strand::Minus)
    /// #     .build()
    /// #     .unwrap()
    /// # );
    /// assert!(transcripts.by_name("NM_001203247.2").is_some());
    /// assert!(transcripts.by_name("invalid_name").is_none());
    /// ```
    pub fn by_name(&self, name: &str) -> Option<&Transcript> {
        match self.name.get(name) {
            Some(id) => self.list.get(*id),
            None => None,
        }
    }

    /// Retrieve all [`Transcript`]s of a gene
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use atg::models::{TranscriptBuilder, Transcripts};
    /// # let mut transcripts = Transcripts::new();
    /// # transcripts.push(TranscriptBuilder::new()
    /// #     .name("NM_001203247.2")
    /// #     .chrom("chr7")
    /// #     .gene("EZH2")
    /// #     .strand(atg::models::Strand::Minus)
    /// #     .build()
    /// #     .unwrap()
    /// # );
    /// assert!(transcripts.by_gene("EZH2").is_some());
    /// assert_eq!(transcripts.by_gene("EZH2").unwrap().len(), 1);
    /// assert!(transcripts.by_gene("Invalid-name").is_none());
    /// ```
    pub fn by_gene(&self, gene: &str) -> Option<Vec<&Transcript>> {
        match self.gene.get(gene) {
            Some(ids) => {
                let mut res:Vec<&Transcript> = Vec::with_capacity(ids.len());
                for id in ids {
                    res.push(&self.list[*id]);
                }
                Some(res)
            },
            None => None,
        }
    }

    /// Add another [`Transcript`]
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::{TranscriptBuilder, Transcripts};
    ///
    /// let mut transcripts = Transcripts::new();
    ///
    /// transcripts.push(TranscriptBuilder::new()
    ///     .name("NM_001203247.2")
    ///     .chrom("chr7")
    ///     .gene("EZH2")
    ///     .strand(atg::models::Strand::Minus)
    ///     .build()
    ///     .unwrap()
    /// );
    ///
    /// transcripts.push(TranscriptBuilder::new()
    ///     .name("NM_001203247.3")
    ///     .chrom("chr7")
    ///     .gene("EZH2")
    ///     .strand(atg::models::Strand::Minus)
    ///     .build()
    ///     .unwrap()
    /// );
    /// assert_eq!(transcripts.len(), 2);
    ///
    /// assert!(transcripts.by_name("NM_001203247.2").is_some());
    /// assert!(transcripts.by_name("NM_001203247.3").is_some());
    /// assert_eq!(transcripts.by_gene("EZH2").unwrap().len(), 2);
    /// ```
    pub fn push(&mut self, record: Transcript) {
        if self.name.get(record.name()).is_some() {
            panic!("A transcript with the name {} already exists", record.name())
        }
        let idx = self.list.len();
        self.name.insert(record.name().to_string(), idx);
        match self.gene.get_mut(record.gene()) {
            Some(x) => x.push(idx),
            None => {
                self.gene.insert(record.gene().to_string(), vec![idx]);
            }
        }
        self.list.push(record);
    }

    /// Returns the number of [`Transcript`]s in the object
    pub fn len(&self) -> usize {
        self.list.len()
    }

    /// Returns true if the object contains no transcripts.
    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    /// Returns a vector of [`Transcript`]s
    pub fn as_vec(&self) -> &Vec<Transcript> {
        &self.list
    }

    /// Consumes and returns a vector of [`Transcript`]s
    pub fn to_vec(self) -> Vec<Transcript> {
        self.list
    }
}

impl Default for Transcripts {
    fn default() -> Self {
        Self::new()
    }
}

impl IntoIterator for Transcripts {
    type Item = Transcript;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.list.into_iter()
    }
}
