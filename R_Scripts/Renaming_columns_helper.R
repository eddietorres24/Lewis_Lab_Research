# Rename a column across *everything* in your environment
rename_everywhere <- function(
    from = "V10",
    to = "gene_id",
    env = .GlobalEnv,
    recurse_lists = TRUE,
    include = NULL,    # regex: only touch objects with names matching this
    exclude = NULL,    # regex: skip objects with names matching this
    verbose = TRUE,
    dry_run = FALSE
){
  # --- helpers ---
  is_s4 <- function(x) methods::is(x, "S4")
  
  fix_df <- function(df) {
    if (inherits(df, "data.table")) {
      if (from %in% names(df)) data.table::setnames(df, from, to)
      return(df)
    }
    if (is.data.frame(df) && from %in% names(df)) {
      names(df)[names(df) == from] <- to
    }
    df
  }
  
  fix_gr <- function(gr) {
    if (!methods::is(gr, "GRanges")) return(gr)
    mc <- S4Vectors::mcols(gr)
    if (!is.null(mc) && ncol(mc) > 0 && from %in% colnames(mc)) {
      cn <- colnames(mc)
      colnames(mc)[cn == from] <- to
      S4Vectors::mcols(gr) <- mc
    }
    gr
  }
  
  fix_any <- function(x) {
    # data.frame/tibble/data.table first (they are lists too)
    if (is.data.frame(x)) return(fix_df(x))
    # GRanges (S4)
    if (methods::is(x, "GRanges")) return(fix_gr(x))
    # Skip environments and many S4 containers
    if (is.environment(x)) return(x)
    if (is_s4(x)) return(x)
    
    # Recurse into plain lists only
    if (is.list(x) && recurse_lists) {
      # Use safe lapply; don't mutate in place
      out <- x
      idx <- seq_along(x)
      for (i in idx) {
        out[[i]] <- tryCatch(
          fix_any(x[[i]]),
          error = function(e) {
            if (verbose) message("  [skip nested] ", conditionMessage(e))
            x[[i]]
          }
        )
      }
      return(out)
    }
    
    x
  }
  
  changed <- character(0)
  objs <- ls(envir = env, all.names = TRUE)
  for (nm in objs) {
    if (!is.null(include) && !grepl(include, nm)) next
    if (!is.null(exclude) &&  grepl(exclude, nm)) next
    
    obj <- get(nm, envir = env)
    new_obj <- tryCatch(fix_any(obj), error = function(e) obj)
    
    # Detect change
    changed_here <- FALSE
    if (!identical(obj, new_obj)) {
      changed_here <- TRUE
      if (!dry_run) assign(nm, new_obj, envir = env)
      changed <- c(changed, nm)
      if (verbose) message(if (dry_run) "[dry] " else "", "Renamed ", from, " -> ", to, " in: ", nm)
    }
    
    # Also try renaming inside lists that contain data.frames but didn’t trigger identical() change
    # (rare, but harmless to skip)
  }
  invisible(changed)
}

###USE HELPER###
rename_everywhere("V10", "gene_id")

