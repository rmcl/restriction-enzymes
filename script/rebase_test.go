package script

import (
	"bufio"
	"strings"
	"testing"

	"github.com/rmcl/restriction-enzymes/enzyme"
)

var bairochSampleInput = `
CC   S=Sigma Chemical Corporation (5/21)
CC   V=Vivantis Technologies (1/18)
CC   X=EURx Ltd. (5/23)
CC   Y=SinaClon BioScience Co. (1/18)
CC
CC

ID   AaaI
ET   R2
AC   RB00001;
OS   Acetobacter aceti ss aceti
PT   XmaIII
RS   CGGCCG, 1;
CR   .
RN   [1]
RA   Tagami H., Tayama K., Tohyama T., Fukaya M., Okumura H., Kawamura Y.,
RA   Horinouchi S., Beppu T.;
RL   FEMS Microbiol. Lett. 56:161-166(1988).
//
`

func TestGetNextBairochRecord(t *testing.T) {
	// Create a String IO with the following content
	input := bairochSampleInput

	record, err := getNextBairochRecord(bufio.NewScanner(strings.NewReader(input)))
	if err != nil {
		t.Fatalf("Error getting next bairoch record: %v", err)
	}

	if recordCC, ok := record["CC"]; ok {
		t.Fatalf("Expected CC to be empty, got %s", recordCC)
	}

	if record["ID"][0] != "AaaI" {
		t.Fatalf("Expected ID to be AaaI, got %s", record["ID"])
	}

	if record["ET"][0] != "R2" {
		t.Fatalf("Expected ET to be R2, got %s", record["ET"])
	}

	if record["AC"][0] != "RB00001;" {
		t.Fatalf("Expected AC to be RB00001;, got %s", record["AC"])
	}

	if record["RA"][0] != "Tagami H., Tayama K., Tohyama T., Fukaya M., Okumura H., Kawamura Y.," {
		t.Fatalf("Expected RA to be Tagami H., Tayama K., Tohyama T., Fukaya M., Okumura H., Kawamura Y.,, got %s", record["RA"])
	}
	if record["RA"][1] != "Horinouchi S., Beppu T.;" {
		t.Fatalf("Expected RA to be Horinouchi S., Beppu T.;, got %s", record["RA"])
	}

}

func TestProcessEnzymeFile(t *testing.T) {
	// Create a String IO with the following content
	input := `
	EcoRI   GAATTC  6       2       0       1       5       0       0
	AanI    TTATAA  6       2       1       3       3       0       0
	AarI    CACCTGC 7       2       0       11      15      0       0
	AasI    GACNNNNNNGTC    12      2       0       7       5       0       0
	BamHI   GGATCC  6       2       0       1       5       0       0
	`

	enzymes := make(map[string]enzyme.Enzyme)
	err := processEnzymeFile(strings.NewReader(input), &enzymes)
	if err != nil {
		t.Fatalf("Error processing enzyme file: %v", err)
	}

	if len(enzymes) != 5 {
		t.Fatalf("Expected 5 enzymes, got %d", len(enzymes))
	}

	if EcoRI, ok := enzymes["EcoRI"]; ok {
		if EcoRI.Name != "EcoRI" {
			t.Fatalf("Expected EcoRI enzyme to have name EcoRI, got %s", EcoRI.Name)
		}
		if EcoRI.Site != "GAATTC" {
			t.Fatalf("Expected EcoRI enzyme to have recognition site GAATTC, got %s", EcoRI.Site)
		}
		if EcoRI.Length != 6 {
			t.Fatalf("Expected EcoRI enzyme to have length 6, got %d", EcoRI.Length)
		}
		if EcoRI.NumberOfCuts != 2 {
			t.Fatalf("Expected EcoRI enzyme to have 2 cuts, got %d", EcoRI.NumberOfCuts)
		}

		if EcoRI.FivePrimeCutSite != 1 {
			t.Fatalf("Expected EcoRI enzyme to have 5' cut site at 1, got %d", EcoRI.FivePrimeCutSite)
		}
		if EcoRI.ThreePrimeCutSite != 5 {
			t.Fatalf("Expected EcoRI enzyme to have 3' cut site at 5, got %d", EcoRI.ThreePrimeCutSite)
		}

		if EcoRI.CutType != enzyme.StickyEnd {
			t.Fatalf("Expected EcoRI enzyme to have sticky end, got %s", EcoRI.CutType)
		}

	} else {
		t.Fatalf("Expected EcoRI enzyme to be present")
	}
}
