import typing as t

from pydantic import BaseModel, Field


class CellTypeSummaryStatisticsResults(BaseModel):
    """
    Represents the data object returned by the CAS API for nearest neighbor annotations.
    """

    class DatasetStatistics(BaseModel):
        dataset_id: str = Field(
            description="The ID of the dataset containing cells", examples=["a7a92fb49-50741b00a-244955d47"]
        )
        count_per_dataset: int = Field(description="The number of cells found in the dataset", examples=[10])
        min_distance: float = Field(
            description="The minimum distance between the query cell and the dataset cells",
            examples=[1589.847900390625],
        )
        max_distance: float = Field(
            description="The maximum distance between the query cell and the dataset cells",
            examples=[1840.047119140625],
        )
        median_distance: float = Field(
            description="The median distance between the query cell and the dataset cells", examples=[1791.372802734375]
        )
        mean_distance: float = Field(
            description="The mean distance between the query cell and the dataset cells", examples=[1791.372802734375]
        )

    class SummaryStatistics(BaseModel):
        cell_type: str = Field(description="The cell type of the cluster of cells", examples=["erythrocyte"])
        cell_count: int = Field(description="The number of cells in the cluster", examples=[94])
        min_distance: float = Field(
            description="The minimum distance between the query cell and the cluster cells",
            examples=[1589.847900390625],
        )
        p25_distance: float = Field(
            description="The 25th percentile distance between the query cell and the cluster cells",
            examples=[1664.875244140625],
        )
        median_distance: float = Field(
            description="The median distance between the query cell and the cluster cells", examples=[1791.372802734375]
        )
        p75_distance: float = Field(
            description="The 75th percentile distance between the query cell and the cluster cells",
            examples=[1801.3585205078125],
        )
        max_distance: float = Field(
            description="The maximum distance between the query cell and the cluster cells",
            examples=[1840.047119140625],
        )
        dataset_ids_with_counts: t.Optional[t.List["CellTypeSummaryStatisticsResults.DatasetStatistics"]] = None

    class NeighborhoodAnnotation(BaseModel):
        """
        Represents the data object returned by the CAS API for a single nearest neighbor annotation.
        """

        query_cell_id: str = Field(description="The ID of the querying cell", examples=["ATTACTTATTTAGTT-12311"])
        matches: t.List["CellTypeSummaryStatisticsResults.SummaryStatistics"]

    data: t.List["CellTypeSummaryStatisticsResults.NeighborhoodAnnotation"] = Field(description="The annotations found")


CellTypeSummaryStatisticsResults.model_rebuild()


class CellTypeOntologyAwareResults(BaseModel):
    """
    Represents the data object returned by the CAS API for a ontology-aware annotations.
    """

    class Match(BaseModel):
        score: float = Field(description="The score of the match", examples=[0.789])
        cell_type_ontology_term_id: str = Field(
            description="The ontology term ID of the cell type for the match", examples=["CL:0000121"]
        )
        cell_type: str = Field(description="The cell type of the match", examples=["erythrocyte"])

    class OntologyAwareAnnotation(BaseModel):
        """
        Represents the data object returned by the CAS API for a single ontology-aware annotation.
        """

        query_cell_id: str = Field(description="The ID of the querying cell", examples=["ATTACTTATTTAGTT-12311"])
        matches: t.List["CellTypeOntologyAwareResults.Match"] = Field(
            description="The matches found for the querying cell"
        )
        total_weight: float = Field(description="The total weight of the matches", examples=[11.23232])
        total_neighbors: int = Field(description="The total number of neighbors matched", examples=[1023])
        total_neighbors_unrecognized: int = Field(
            description="The total number of neighbors that were not recognized", examples=[5]
        )

    data: t.List["CellTypeOntologyAwareResults.OntologyAwareAnnotation"] = Field(description="The annotations found")


CellTypeOntologyAwareResults.model_rebuild()


class MatrixQueryResults(BaseModel):
    """
    Represents the data object returned by the CAS API when performing a cell matrix query
    (e.g. a query of the cell database using a matrix).
    """

    class Match(BaseModel):
        cas_cell_index: float = Field(description="CAS-specific ID of a single cell", examples=[123])
        distance: float = Field(
            description="The distance between this querying cell and the found cell", examples=[0.123]
        )

    class MatrixQueryResult(BaseModel):
        """
        Represents the data object returned by the CAS API for a single cell query.
        """

        query_cell_id: str = Field(description="The ID of the querying cell", examples=["ATTACTTATTTAGTT-12311"])
        neighbors: t.List["MatrixQueryResults.Match"]

    data: t.List["MatrixQueryResults.MatrixQueryResult"] = Field(description="The results of the query")


MatrixQueryResults.model_rebuild()


class CellQueryResults(BaseModel):
    """
    Represents the data object returned by the CAS API for a cell query.
    """

    class CellariumCellMetadata(BaseModel):
        cas_cell_index: int = Field(description="The CAS-specific ID of the cell", examples=[123])
        cell_type: t.Optional[str] = Field(description="The cell type of the cell", examples=["enterocyte"])
        assay: t.Optional[str] = Field(description="The assay used to generate the cell", examples=["10x 3' v2"])
        disease: t.Optional[str] = Field(description="The disease state of the cell", examples=["glioblastoma"])
        donor_id: t.Optional[str] = Field(description="The ID of the donor of the cell", examples=["H20.33.013"])
        is_primary_data: t.Optional[bool] = Field(description="Whether the cell is primary data", examples=[True])
        development_stage: t.Optional[str] = Field(
            description="The development stage of the cell donor", examples=["human adult stage"]
        )
        organism: t.Optional[str] = Field(description="The organism of the cell", examples=["Homo sapiens"])
        self_reported_ethnicity: t.Optional[str] = Field(
            description="The self reported ethnicity of the cell donor", examples=["Japanese"]
        )
        sex: t.Optional[str] = Field(description="The sex of the  cell donor", examples=["male"])
        suspension_type: t.Optional[str] = Field(description="The cell suspension types used", examples=["nucleus"])
        tissue: t.Optional[str] = Field(
            description="The tissue-type that the cell was a part of", examples=["cerebellum"]
        )
        total_mrna_umis: t.Optional[int] = Field(
            description="The count of mRNA UMIs associated with this cell", examples=[24312]
        )

        # Ontology term IDs for the fields
        cell_type_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the type of the cell", examples=["CL:0000121"]
        )
        assay_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the assay used to generate the cell", examples=["EFO:0010550"]
        )
        disease_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the disease state of the cell", examples=["PATO:0000461"]
        )
        development_stage_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the development stage of the cell donor",
            examples=["HsapDv:0000053"],
        )
        organism_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the organism of the cell", examples=["NCBITaxon:9606"]
        )
        self_reported_ethnicity_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the self reported ethnicity of the cell donor",
            examples=["HANCESTRO:0019"],
        )
        sex_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the sex of the cell donor", examples=["PATO:0000384"]
        )
        tissue_ontology_term_id: t.Optional[str] = Field(
            description="The ID used by the ontology for the tissue type that the cell was a part of",
            examples=["UBERON:0002037"],
        )

    data: t.List["CellQueryResults.CellariumCellMetadata"] = Field(description="The metadata of the found cells")


CellQueryResults.model_rebuild()
