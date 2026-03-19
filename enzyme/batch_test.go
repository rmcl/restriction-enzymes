package enzyme

import (
	"regexp"
	"testing"

	"github.com/rmcl/restriction-enzymes/constants"
)

func TestGetNextRecognitionSiteForBatch(unittest *testing.T) {

	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		FIXTURES["EcoRI"],
	)
	sequence := "AATAGACAGAATTCGATTCACCAGAGGTCTCATAGACAAACCAGAGAAAAAAAA"

	results := batch.GetNextRecognitionSite(sequence, 0, false)
	if (results[0].RecognitionSiteIndex != 8) || results[0].Enzyme.Name != "EcoRI" || (results[0].Strand != constants.Watson) {
		unittest.Errorf("Expected 8, nil, EcoRI, got %d, %v, %s", results[0].RecognitionSiteIndex, results[0].Enzyme, results[0].Strand)
	}

	results = batch.GetNextRecognitionSite(sequence, 9, false)
	if (results[0].RecognitionSiteIndex != 25) || results[0].Enzyme.Name != "BsaI" || (results[0].Strand != constants.Watson) {
		unittest.Errorf("Expected 25, nil, BsaI, got %d, %v, %s", results[0].RecognitionSiteIndex, results[0].Enzyme, results[0].Strand)
	}
}

func TestAddEnzymeToBatch(unittest *testing.T) {

	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		FIXTURES["EcoRI"],
	)

	if len(batch.Enzymes) != 2 {
		unittest.Errorf("Expected 2 enzymes, got %d", len(batch.Enzymes))
	}

	batch.Add(FIXTURES["BamHI"])

	if len(batch.Enzymes) != 3 {
		unittest.Errorf("Expected 3 enzymes, got %d", len(batch.Enzymes))
	}

}

func TestSearchBatch(unittest *testing.T) {
	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		FIXTURES["EcoRI"],
	)
	sequence := "AATAGACAGAATTCGATTCACCAGAGGTCTCATAGACAAACCAGAGAAAAAAAA"

	results, err := batch.Search(sequence, false)
	if err != nil {
		unittest.Errorf("Error searching: %v", err)
	}

	if len(results) != 2 {
		unittest.Errorf("Expected 2 enzymes, got %d", len(results))
	}

	if len(results["EcoRI"]) != 1 {
		unittest.Errorf("Expected 1 EcoRI site, got %d", len(results["EcoRI"]))
	}

	if len(results["BsaI"]) != 1 {
		unittest.Errorf("Expected 1 BsaI site, got %d", len(results["BsaI"]))
	}

}

func TestSearchBatchPrefersLongerOverlappingRecognitionSite(unittest *testing.T) {
	shortOverlap := Enzyme{
		Name:              "ShortOverlap",
		Site:              "GTCTC",
		Length:            5,
		Substrate:         "DNA",
		RegexpFor:         regexp.MustCompile("(?i)GTCTC"),
		RegexpRev:         regexp.MustCompile("(?i)GAGAC"),
		NumberOfCuts:      2,
		CutType:           "sticky",
		FivePrimeCutSite:  6,
		ThreePrimeCutSite: 10,
	}

	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		shortOverlap,
	)

	sequence := "AAAGGTCTCAAAGAGACCTTT"

	results, err := batch.Search(sequence, false)
	if err != nil {
		unittest.Fatalf("Error searching: %v", err)
	}

	if len(results["BsaI"]) != 2 {
		unittest.Fatalf("Expected 2 BsaI sites, got %d (%v)", len(results["BsaI"]), results["BsaI"])
	}

	if len(results["ShortOverlap"]) != 2 {
		unittest.Fatalf("Expected 2 ShortOverlap sites, got %d (%v)", len(results["ShortOverlap"]), results["ShortOverlap"])
	}
}

func TestGetNextRecognitionSiteReturnsAllOverlappingMatchesAtPosition(unittest *testing.T) {
	shortOverlap := Enzyme{
		Name:              "ShortOverlap",
		Site:              "GTCTC",
		Length:            5,
		Substrate:         "DNA",
		RegexpFor:         regexp.MustCompile("(?i)GTCTC"),
		RegexpRev:         regexp.MustCompile("(?i)GAGAC"),
		NumberOfCuts:      2,
		CutType:           "sticky",
		FivePrimeCutSite:  6,
		ThreePrimeCutSite: 10,
	}

	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		shortOverlap,
	)

	results := batch.GetNextRecognitionSite("AAAGAGACCTTT", 0, false)
	if len(results) != 2 {
		unittest.Fatalf("Expected 2 overlapping matches, got %d", len(results))
	}

	if results[0].Enzyme.Name != "BsaI" || results[0].Strand != constants.Crick {
		unittest.Fatalf("Expected first overlapping match to be BsaI on the crick strand, got %s %s", results[0].Enzyme.Name, results[0].Strand)
	}

	if results[1].Enzyme.Name != "ShortOverlap" || results[1].Strand != constants.Crick {
		unittest.Fatalf("Expected second overlapping match to be ShortOverlap on the crick strand, got %s %s", results[1].Enzyme.Name, results[1].Strand)
	}
}
