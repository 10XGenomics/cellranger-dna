// Package vttime contains the definitions and implementations for the Vitess
// time library. This package is based on Google's TrueTime, as described
// in this Spanner paper for instance:
// http://static.googleusercontent.com/media/research.google.com/en//archive/spanner-osdi2012.pdf
//
// The idea is that a timestamp is not enough, as clocks will drift
// apart between computers. However, it is usually possible to know
// how much drift happens. So instead of returning a timestamp that
// may be wrong, we return an interval [earliest, latest] with the
// following guarantees:
// - current time is greater or equal to 'earliest'.
// - current time is less or equal to 'latest'.
//
// When comparing two intervals, we know one of them is smaller if
// there is no overlap and it is before:
//   [--------]
//                 [-----------]
// If there is overlap, we can't say for sure:
//   [--------]
//        [----------]
//
// However, if the goal is to be sure we are producing events that
// clients will know are chonologically ordered, it is then possible
// to sleep for a few milliseconds and guarantee that. This becomes
// handy in Paxos-like algorithms, for instance. See the paper for
// more details.
package vttime
