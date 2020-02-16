//! # Supermer
//! supermer provides a simple method for calculating
//! supermers from a sequence of hash values in a streaming fashion.
//!
//! In this routine a supermer is defined as a substring with length
//! of at least `w`, containing a kmer which is the smaller than all other kmers
//! in the substring. The supermer may be a extending up to `w` either size
//! of the minimal kmer.
//!
//! ## Graphically
//!
//! ```text
//!                  min       Head      End
//!                  <---------(<=w)------->
//!                  <--(k)-->
//!   ===============|=======|============== ---> Growth
//!   <---------(<=w)-------->
//!   Start     Tail
//! ```
//!
//! ## Rules
//! The tail and head must each never exceed `w` in length.
//! The length of the supermer is between `w` and 2`w` - `k`.
//!
//! ## Note
//! The invariant that length must be at least `w` can only be
//! maintained by defining the right most `kmer` as the smallest
//! in the case of a tie.

use std::cmp::Ordering;
/// # Supermer
///
/// Thus the overall length may never exceed 2w - k.
/// The overall length must always be at least w.
#[derive(Debug, PartialEq)]
pub struct Supermer<H> {
    start: usize,
    end: usize,
    min: Minima<H>,
}

/// A Minima is a combination of position and value.
#[derive(Debug, Clone, PartialEq, Copy)]
pub struct Minima<H> {
    pos: usize,
    min: H,
}

/// Minima ordering is by value (and technically position,
/// in case of a draw).
impl<H: PartialOrd + Copy> std::cmp::PartialOrd for Minima<H> {
    fn partial_cmp(&self, rhs: &Minima<H>) -> Option<Ordering> {
        self.min.partial_cmp(&rhs.min)
    }
}

/// # Tracker
///
/// In order to efficiently create supermers in a
/// streaming fashion, we must maintain a list of
/// three minima. The first is the proposed minima
/// of the current supermer, the next if present, is the proposed minima
/// of the next supermer. The final minima, if present
/// is the proposed minima to the successor of the next supermer.
#[derive(PartialEq, Debug, Default)]
struct Tracker<H> {
    // Better variable names?
    a: Option<Minima<H>>,
    b: Option<Minima<H>>,
    c: Option<Minima<H>>,
}

impl<H: PartialOrd + Copy + Default> Tracker<H> {
    fn track(&mut self, val: Minima<H>) {
        if let Some(z) = &mut self.a {
            if val <= *z {
                *z = val;
                self.b = None;
                self.c = None;
                return;
            }
        } else {
            self.a = Some(val);
            return;
        }

        if let Some(z) = &mut self.b {
            if val <= *z {
                *z = val;
                self.c = None;
                return;
            }
        } else {
            self.b = Some(val);
            return;
        }

        if let Some(z) = &mut self.c {
            if val <= *z {
                *z = val;
            }
        } else {
            self.c = Some(val)
        }
    }

    fn is_minima(&self, val: Minima<H>) -> bool {
        if let Some(m) = self.a {
            val <= m
        } else {
            true
        }
    }

    fn minima(&self) -> Minima<H> {
        self.a.expect("Expected a minima!")
    }

    fn reset(&mut self) {
        self.a = self.b;
        self.b = self.c;
        self.c = None;
        assert!(self.a.is_some())
    }
}

#[derive(Debug, PartialEq)]
pub struct SupermerBuilder<H> {
    k: usize,
    w: usize,
    start: usize,
    end: usize,
    partial: bool,
    tracker: Tracker<H>,
}

impl<H: PartialOrd + Copy + Default> SupermerBuilder<H> {
    pub fn new(k: usize, w: usize) -> Self {
        assert!(k <= w);
        Self {
            k,
            w,
            start: 0,
            end: 0,
            partial: true,
            tracker: Tracker::default(),
        }
    }

    fn min_end(&self) -> usize {
        self.end + self.k
    }

    fn in_head(&self) -> bool {
        if self.end != 0 {
            self.min_end() - self.tracker.minima().pos < self.w
        } else {
            true
        }
    }

    fn in_tail(&self) -> bool {
        self.min_end() - self.start <= self.w
    }

    fn res(&self) -> Supermer<H> {
        Supermer {
            start: self.start,
            end: self.min_end(),
            min: self.tracker.minima(),
        }
    }

    fn track(&mut self, loc: Minima<H>) {
        self.tracker.track(loc);
        self.end += 1;
    }

    /// Given a value, apply it to the under-construction
    /// supermer. Returns a supermer if one was completed.
    pub fn build(&mut self, val: H) -> Option<Supermer<H>> {
        let loci = Minima {
            min: val,
            pos: self.end,
        };
        if self.tracker.is_minima(loci) && !self.in_tail() {
            let ret = self.res();
            self.start = self.min_end() - self.w;
            self.track(loci);
            self.partial = false;
            return Some(ret);
        } else if !self.in_head() {
            let ret = self.res();
            self.start = self.tracker.minima().pos + 1;
            self.tracker.reset();
            self.track(loci);
            self.partial = false;
            return Some(ret);
        }

        self.track(loci);
        self.partial = true;
        None
    }

    /// Return any supermer if there is one unfinished.
    pub fn finish(&mut self) -> Option<Supermer<H>> {
        if self.partial {
	    self.partial = false;
            Some(self.res())
        } else {
            None
        }
    }
}

/// Given an iterator which produces a stream of kmer based values `H`,
/// computes a stream of Supermers.
pub struct SupermerIter<H, I> {
    iter: I,
    builder: SupermerBuilder<H>
}

impl<H: PartialOrd + Copy + Default, I: Iterator<Item=H>> SupermerIter<H, I> {
    /// Create a new SupermerIter from a kmer iterator. The argument `k`
    /// should match that of the kmer iterator.
    pub fn new(iter: I, k: usize, w: usize) -> Self {
	Self{ iter, builder: SupermerBuilder::new(k, w) }
    }
}

impl<H: PartialOrd + Copy + Default, I: Iterator<Item=H>> Iterator for SupermerIter<H, I> {
    type Item = Supermer<H>;
    fn next(&mut self) -> Option<Supermer<H>> {
	while let Some(h) = self.iter.next() {
	    if let Some(s) = self.builder.build(h) {
		return Some(s)
	    }
	}
	self.builder.finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_tracker() {
        let mut t: Tracker<usize> = Tracker::default();
        assert_eq!(
            t,
            Tracker {
                a: None,
                b: None,
                c: None
            }
        );
        t.track(Minima { pos: 0, min: 10 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 0 }),
                b: None,
                c: None
            }
        );
        t.track(Minima { pos: 1, min: 10 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 1 }),
                b: None,
                c: None
            }
        );

        t.track(Minima { pos: 2, min: 13 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 1 }),
                b: Some(Minima { min: 13, pos: 2 }),
                c: None
            }
        );

        t.track(Minima { pos: 3, min: 12 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 1 }),
                b: Some(Minima { min: 12, pos: 3 }),
                c: None
            }
        );
        t.track(Minima { pos: 4, min: 16 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 1 }),
                b: Some(Minima { min: 12, pos: 3 }),
                c: Some(Minima { min: 16, pos: 4 })
            }
        );

        t.track(Minima { pos: 5, min: 11 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 1 }),
                b: Some(Minima { min: 11, pos: 5 }),
                c: None
            }
        );
        t.track(Minima { pos: 6, min: 14 });
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 10, pos: 1 }),
                b: Some(Minima { min: 11, pos: 5 }),
                c: Some(Minima { min: 14, pos: 6 }),
            }
        );
        t.reset();
        assert_eq!(
            t,
            Tracker {
                a: Some(Minima { min: 11, pos: 5 }),
                b: Some(Minima { min: 14, pos: 6 }),
                c: None
            }
        );
        assert!(t.is_minima(Minima { min: 8, pos: 7 }));
        assert!(!t.is_minima(Minima { min: 18, pos: 7 }));
    }

    // Currently unused...
    fn assert_bounds<H: Copy + PartialOrd>(k: usize, s: usize, sm: &Supermer<H>) {
	assert!(sm.start < sm.end);
	assert!(sm.end - sm.min.pos <= s, "Violated max head length");
	assert!(sm.end - sm.start <= s + s - k, "Violated max length");
	assert!(sm.end - sm.start >= s, "Violated min length");
	assert!(sm.min.pos - sm.start <= s - k, "Violoated max tail length");
    }
}
