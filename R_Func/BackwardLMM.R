#' Retrieve Model Summary, ANOVA, and R-squared via Backward selection
#'
#' @description
#' \code{BackwardLMM} fits a full linear mixed-effects model using \code{lmer} from the \pkg{lme4} package and
#' performs backward model selection by iteratively removing non-significant fixed-effect terms.
#' At each step, the function prioritizes the removal of the highest-order interaction terms first, then considers
#' main effects. When \code{keep_hierarchy = TRUE}, main effects involved in retained interactions are protected
#' (hierarchical principle / marginality).
#'
#' The function returns sample sizes for retained terms, the final model,
#' its summary, a Type III ANOVA table (via \pkg{car}), and the marginal and conditional R-squared values
#' (via \pkg{MuMIn}).
#'
#' @param Data A data frame containing the dataset to be modeled.
#' @param expVar A character string specifying the response variable.
#' @param respExpr A character string specifying the fixed-effects portion of the model (including main effects and interaction terms).
#' @param random A character string specifying the random effects structure (e.g., \code{"(1 | Session/Bac) + (1 | Pop)"}).
#' @param alpha Numeric value giving the significance threshold used to remove terms (default \code{0.05}).
#' @param keep_hierarchy Logical; if \code{TRUE}, main effects that are part of retained interaction terms are not removed.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{SampleSize}{A data frame reporting sample sizes for each retained fixed-effect term.}
#'   \item{ModL}{The final refined mixed effects model (an object of class \code{lmerMod}).}
#'   \item{ModSum}{The summary of the final model.}
#'   \item{ModAnov}{The Type III ANOVA table for the final model.}
#'   \item{Rsquared}{A numeric vector containing the marginal and conditional R-squared values as computed by \code{MuMIn::r.squaredGLMM}.}
#' }
#'
#' @details
#' \code{BackwardLMM} first filters the input data to exclude missing values in the response variable.
#' A full model is then fitted using the formula constructed from the response variable, the specified fixed effects,
#' and the random-effects structure.
#'
#' Backward selection proceeds iteratively:
#' \enumerate{
#'   \item Interaction terms are tested and removed first (highest-order interactions prioritized).
#'   \item Main effects are then tested and removed, unless protected by \code{keep_hierarchy = TRUE}.
#' }
#'
#' The function returns the refined model along with
#' its summary, Type III ANOVA table, and R-squared values.
#'
#' **Dependencies:**
#' This function depends on:
#' \itemize{
#'   \item \strong{lme4} (model fitting with \code{lmer}),
#'   \item \strong{car} (Type III ANOVA with \code{car::Anova}),
#'   \item \strong{MuMIn} (R-squared computation with \code{MuMIn::r.squaredGLMM})
#' }
#'
#' @importFrom lme4 lmer
#' @importFrom car Anova
#' @importFrom MuMIn r.squaredGLMM
#'
#' @author Quentin PETITJEAN [quentin.petitjean@inrae.fr]
#'
#' @date 22/01/2026
#'
#' @export
 
BackwardLMM <- function(Data = NULL, # a dataframe containing the data to model
                       expVar = NULL, # the explainatory variable
                       respExpr = NULL, # the full set of response variables, including interactions terms
                       random = NULL, # the random effect 
                       alpha = 0.05,
                       keep_hierarchy = TRUE)      # keep main effects if they appear in retained interactions 
  {
  
  toMod <- Data[!is.na(Data[[expVar]]), ]
  
  # Build a full model and refined it using AICc
  ModCur <-
    lme4::lmer(
      as.formula(paste0(
        paste0(expVar, " ~ "), respExpr, paste0("+", random)
      )),
      data = toMod,
      na.action = "na.fail",
      REML = T
    )
  
  # helper: number of ":" gives interaction order (":" -> 2-way, "::" -> 3-way, etc.)
  colon_count <- function(x) {
    m <- gregexpr(":", x, fixed = TRUE)[[1]]
    if (m[1] == -1) 0 else length(m)
  }
  
  
  # helper: extract p-values for each term currently in the fixed part
  get_term_pvals <- function(mod) {
    fixed_terms <- attr(stats::terms(lme4::nobars(stats::formula(mod))), "term.labels")
    
    if (length(fixed_terms) == 0) return(setNames(numeric(0), character(0)))
      a <- car::Anova(mod, type = 3)
      pcol <- grep("^Pr\\(", colnames(a), value = TRUE)[1]
      p <- setNames(a[, pcol], rownames(a))
      p <- p[names(p) != "(Intercept)"]
    
    # keep only terms that are truly in the fixed part right now (and keep ordering stable)
    p[intersect(fixed_terms, names(p))]
  }
  
  
  # Backward elimination loop
  repeat {
    fixed_terms <- attr(stats::terms(lme4::nobars(stats::formula(ModCur))), "term.labels")
    if (length(fixed_terms) == 0) break
    
    ord <- 1 + vapply(fixed_terms, colon_count, numeric(1))
    names(ord) <- fixed_terms
    
    pvals <- get_term_pvals(ModCur)
    if (length(pvals) == 0) break
    
    # --- Phase 1: try to drop non-significant interactions (highest order first) ---
    int_terms <- fixed_terms[ord > 1]
    if (length(int_terms) > 0) {
      p_int <- pvals[int_terms]
      p_int <- p_int[!is.na(p_int)]
      
      if (length(p_int) > 0 && any(p_int > alpha)) {
        cand <- names(p_int)[p_int > alpha]
        cand_ord <- ord[cand]
        max_ord <- max(cand_ord)
        cand2 <- cand[cand_ord == max_ord]
        drop_term <- cand2[which.max(p_int[cand2])]
        
        new_form <- stats::update(stats::formula(ModCur), paste(". ~ . -", drop_term))
        ModNew <- tryCatch(
          lme4::lmer(new_form, data = toMod, na.action = "na.fail", REML = FALSE),
          error = function(e) NULL
        )
        if (is.null(ModNew)){ break
        message("Dropping interaction: ", drop_term, " (p=", signif(p_int[drop_term], 3), ")")
        }
        ModCur <- ModNew
        next
      }
    }
    
    # --- Phase 2: drop non-significant main effects not needed for hierarchy ---
    protected <- character(0)
    if (keep_hierarchy && length(int_terms) > 0) {
      protected <- unique(unlist(strsplit(int_terms, ":", fixed = TRUE)))
    }
    
    main_terms <- setdiff(fixed_terms[ord == 1], protected)
    if (length(main_terms) > 0) {
      p_main <- pvals[main_terms]
      p_main <- p_main[!is.na(p_main)]
      
      if (length(p_main) > 0 && any(p_main > alpha)) {
        drop_term <- names(p_main)[which.max(p_main)]
        
        new_form <- stats::update(stats::formula(ModCur), paste(". ~ . -", drop_term))
        ModNew <- tryCatch(
          lme4::lmer(new_form, data = toMod, na.action = "na.fail", REML = FALSE),
          error = function(e) NULL
        )
        if (is.null(ModNew)){ break
          message("Dropping main effect: ", drop_term, " (p=", signif(p_main[drop_term], 3), ")")}
        
        ModCur <- ModNew
        next
      }
    }
    
    # nothing left to drop
    break
  }

  # Final refit with REML=TRUE for inference/reporting
  ModFinal <- lme4::lmer(
    stats::formula(ModCur),
    data = toMod,
    na.action = "na.fail",
    REML = TRUE
  )
  
  # Terms kept (fixed effects only)
  toKeep <- attr(stats::terms(lme4::nobars(stats::formula(ModFinal))), "term.labels")
  
  # retrieve sample size for each variable remaining in the final Model
  sampleSize <- c()
  for(i in seq_along(toKeep)){
    if(toKeep[[i]] %in% names(toMod)){
      keepVar <- toKeep[[i]]
    }else{
      keepVar <-  gsub(":", ".", toKeep[[i]])
    }
    tab <- as.data.frame(table(toMod[[keepVar]]))
    tab[["treat"]] <- rep(keepVar, nrow(tab))
    sampleSize <- rbind(sampleSize, tab) 
  }
  
  
  # retrieve the results in a list
  Res <- list(
    SampleSize = sampleSize,
    ModL = ModFinal,
    ModSum = summary(ModFinal),
    ModAnov = car::Anova(ModFinal, type = "3"),
    Rsquared = MuMIn::r.squaredGLMM(ModFinal)
  )
  
  return(Res)
}
