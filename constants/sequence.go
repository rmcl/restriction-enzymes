package constants

type Strand string

const (
	Watson Strand = "watson"
	Crick  Strand = "crick"
)

type SequenceGeometry string

const (
	Linear   SequenceGeometry = "linear"
	Circular SequenceGeometry = "circular"
)
