#' @title DNA 5-mer one-hot encoding
#' @name dna_encoding
#' @description
#' This function converts DNA sequences into a matrix-like data frame in which
#' each nucleotide position becomes a separate column. Each entry is stored as
#' a factor with four possible levels: `"A"`, `"T"`, `"C"`, `"G"`.
#'
#' This encoding format is often used in machine learning models, where
#' sequence characters must be transformed into numeric or categorical values
#' before prediction.
#'
#' @param dna_strings A character vector of DNA sequences. All sequences must
#' have the same length (e.g. `"ATCGA"` or `c("ATCGA", "TTGCA")`).
#'
#' @return A data frame where:
#' \itemize{
#'   \item Each row corresponds to one DNA sequence
#'   \item Each column corresponds to a nucleotide position
#'   \item Column names follow the format `"nt_pos1"`, `"nt_pos2"`, etc.
#'   \item Values are factors with levels `"A"`, `"T"`, `"C"`, `"G"`
#' }
#'
#' @examples
#' # Example with one sequence
#' dna_encoding("ATCGA")
#'
#' # Example with multiple sequences
#' dna_encoding(c("ATCGA", "TTGCA"))
#'
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' @title Batch prediction for m6A
#' @name prediction_multiple
#' @importFrom stats predict
#' @import randomForest
#' @description
#' This function predicts m6A methylation probability and binary status
#' for a batch of RNA sites using a trained classification model.
#' It validates the input column structure, encodes categorical variables,
#' performs one-hot encoding for 5-mer sequences, and returns prediction results.
#'
#' @param ml_fit A trained classification model that supports
#' `predict(..., type = "prob")`, such as a tidymodels or caret model,
#' with a `"Positive"` probability column.
#'
#' @param feature_df A data frame containing features of RNA sites to predict.
#' Required columns:
#' \itemize{
#'   \item `gc_content` – numeric GC content (0–1)
#'   \item `RNA_type` – one of `"mRNA"`, `"lincRNA"`, `"lncRNA"`, `"pseudogene"`
#'   \item `RNA_region` – one of `"CDS"`, `"intron"`, `"3'UTR"`, `"5'UTR"`
#'   \item `exon_length` – numeric exon length
#'   \item `distance_to_junction` – numeric distance to exon–intron boundary
#'   \item `evolutionary_conservation` – numeric conservation score
#'   \item `DNA_5mer` – a string of 5 nucleotides (e.g., `"ATGCG"`)
#' }
#'
#' @param positive_threshold Probability threshold for assigning `"Positive"` status.
#' Default is `0.5`.
#'
#' @return A data frame containing the original features plus:
#' \itemize{
#'   \item `predicted_m6A_prob` – predicted probability (numeric)
#'   \item `predicted_m6A_status` – `"Positive"` or `"Negative"` (character)
#' }
#'
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' input_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' prediction_multiple(model, input_df, positive_threshold = 0.6)
#'
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  prob <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- prob
  feature_df$predicted_m6A_status <- ifelse(prob > positive_threshold, "Positive", "Negative")
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

#' @title Single-site m6A prediction
#' @name prediction_single
#' @description
#' This function predicts whether a single RNA site is likely to be **m6A-modified**
#' based on its sequence features and genomic context.
#' It creates a one-row data frame using the input arguments and passes it
#' to `prediction_multiple()` to reuse the same encoding and prediction logic.
#'
#' @param ml_fit A trained machine learning model that supports `predict(..., type = "prob")`.
#'
#' @param gc_content Numeric value indicating the GC content of the region.
#'
#' @param RNA_type Character string: one of `"mRNA"`, `"lincRNA"`, `"lncRNA"`, or `"pseudogene"`.
#'
#' @param RNA_region Character string: one of `"CDS"`, `"intron"`, `"3'UTR"`, or `"5'UTR"`.
#'
#' @param exon_length Numeric value indicating the exon length at the site.
#'
#' @param distance_to_junction Numeric value for distance to nearest exon–intron junction.
#'
#' @param evolutionary_conservation Numeric score for evolutionary conservation (e.g., phastCons).
#'
#' @param DNA_5mer A 5-mer DNA sequence (e.g., `"ATGCG"`).
#'
#' @param positive_threshold Probability threshold to call a site "Positive". Default is 0.5.
#'
#' @return A named vector with:
#' \itemize{
#'   \item `predicted_m6A_prob` – Predicted probability of m6A
#'   \item `predicted_m6A_status` – Either `"Positive"` or `"Negative"`
#' }
#'
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(
#'   ml_fit = model,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 120,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATCGAT",
#'   positive_threshold = 0.5
#' )
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5) {
  input_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  result_df <- prediction_multiple(ml_fit, input_df, positive_threshold)
  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob,
    predicted_m6A_status = result_df$predicted_m6A_status
  )
  return(returned_vector)
}
