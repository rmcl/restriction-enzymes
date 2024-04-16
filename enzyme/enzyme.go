package enzyme

import "regexp"

type EnzymeNumberOfCuts int

const (
	OneCut      EnzymeNumberOfCuts = 1
	TwoCuts     EnzymeNumberOfCuts = 2
	UnknownCuts EnzymeNumberOfCuts = 0
)

type EnzymeCutType string

const (
	BluntEnd  EnzymeCutType = "blunt"
	StickyEnd EnzymeCutType = "sticky"
)

type Enzyme struct {
	Name string

	Site      string
	Length    int
	Substrate string

	RegexpFor *regexp.Regexp
	RegexpRev *regexp.Regexp

	CutType EnzymeCutType

	OverhangLength   int
	OverhangSequence string

	NumberOfCuts EnzymeNumberOfCuts

	FivePrimeCutSite  int
	ThreePrimeCutSite int

	FivePrimeCutSite2  int
	ThreePrimeCutSite2 int

	RebaseId                int
	InactivationTemperature int
	OptimalTemperature      int

	Uri string

	Suppliers []string

	References []string
}
