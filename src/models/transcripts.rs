use std::collections::HashMap;

use crate::models::Transcript;

/// A convinience wrapper to handle  large amounts of `Transcript`s
pub struct Transcripts {
    list: Vec<Transcript>,
    name: HashMap<String, usize>,
    gene: HashMap<String, usize>,
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

    pub fn by_name(&self, name: &str) -> Option<&Transcript> {
        match self.name.get(name) {
            Some(id) => self.list.get(*id),
            None => None,
        }
    }

    pub fn by_gene(&self, gene: &str) -> Option<&Transcript> {
        match self.gene.get(gene) {
            Some(id) => self.list.get(*id),
            None => None,
        }
    }

    pub fn push(&mut self, record: Transcript) {
        let idx = self.list.len();
        self.name.insert(record.name().to_string(), idx);
        self.name.insert(record.gene().to_string(), idx);
        self.list.push(record);
    }

    pub fn len(&self) -> usize {
        self.list.len()
    }

    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    pub fn as_vec(&self) -> &Vec<Transcript> {
        &self.list
    }
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
