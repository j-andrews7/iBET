#' Swap Data Frame Row Names Based on ID Mapping
#'
#' This internal utility function replaces the row names of a data frame (`df`)
#' with new identifiers based on a provided mapping (`id_map`). It matches the
#' current row names of `df` to a vector of original IDs (`orig_ids`), retains
#' only the rows with valid matches, and assigns new row names from `id_map`.
#' The original row names are preserved in a new column called `ORIGINAL_ROWS`.
#'
#' @param df A data frame whose row names are to be swapped.
#' @param id_map A named character vector or list mapping old row names to new row names.
#' @param orig_ids A character vector of original IDs to match against the row names of `df`.
#'
#' @return A data frame with updated row names and an additional column
#'   `ORIGINAL_ROWS` containing the original row names. If no matches are found,
#'   returns an empty data frame with the same structure as `df`.
#'
#' @keywords internal
#' @examples
#' # Not intended for direct use.
.swap_res_rownames <- function(df, id_map, orig_ids) {
    # Match result rownames to mat original rownames
    matched_indices <- match(rownames(df), orig_ids)
    valid_matches <- !is.na(matched_indices)

    if (sum(valid_matches) == 0) {
        warning("No matching features found between results and expression matrix after swap.rownames")
        return(df[0, , drop = FALSE]) # Return empty dataframe with same structure
    }

    # Keep only matching rows
    df_matched <- df[valid_matches, , drop = FALSE]

    # Store original rownames before swapping
    df_matched$ORIGINAL_ROWS <- rownames(df_matched)

    # Assign new rownames
    new_rownames <- id_map[rownames(df)[valid_matches]]
    rownames(df_matched) <- new_rownames

    return(df_matched)
}
