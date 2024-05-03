package enzyme

import (
	"fmt"
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
	if (results[0].Position != 8) || results[0].Enzyme.Name != "EcoRI" || (results[0].Strand != constants.Watson) {
		unittest.Errorf("Expected 8, nil, EcoRI, got %d, %v, %s", results[0].Position, results[0].Enzyme, results[0].Strand)
	}

	results = batch.GetNextRecognitionSite(sequence, 9, false)
	if (results[0].Position != 25) || results[0].Enzyme.Name != "BsaI" || (results[0].Strand != constants.Watson) {
		unittest.Errorf("Expected 25, nil, BsaI, got %d, %v, %s", results[0].Position, results[0].Enzyme, results[0].Strand)
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

func TestSearch(unittest *testing.T) {
	batch := NewRestrictionBatch(
		FIXTURES["BsaI"],
		FIXTURES["EcoRI"],
	)
	sequence := "AATAGACAGAATTCGATTCACCAGAGGTCTCATAGACAAACCAGAGAAAAAAAA"

	results, err := batch.Search(sequence, false)
	if err != nil {
		unittest.Errorf("Error searching: %v", err)
	}

	fmt.Println(results)
	unittest.Errorf("Not implemented")
}
