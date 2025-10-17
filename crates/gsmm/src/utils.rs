
// Used to enable using newtypes for indexing, even though I am using Vecs internally
pub trait GsmmKey {
    fn to_index(&self) -> usize;
}

pub struct GsmmMap<K, V> {
    values: Vec<V>,
    _marker: std::marker::PhantomData<K>,
}

impl<K: GsmmKey, V> std::ops::Index<K> for GsmmMap<K, V> {
    type Output = V;

    fn index(&self, index: K) -> &Self::Output {
        &self.values[index.to_index()]
    }
}

impl<K: GsmmKey, V> GsmmMap<K, V> {
    pub fn get(&self, key: K) -> Option<&V> {
        self.values.get(key.to_index())
    }
}

pub struct GeneId {
    pub(crate) ind: usize,
}

pub struct MetaboliteIndex {
    pub(crate) ind: usize,
}

pub struct ReactionIndex {
    pub(crate) ind: usize,
}

impl GsmmKey for GeneId {
    fn to_index(&self) -> usize {
        self.ind
    }
}

impl GsmmKey for MetaboliteIndex {
    fn to_index(&self) -> usize {
        self.ind
    }
}

impl GsmmKey for ReactionIndex {
    fn to_index(&self) -> usize {
        self.ind
    }
}