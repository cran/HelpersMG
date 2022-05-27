.getint <- function (cf, parm, ses, level = 0.95, df = NULL, ...) 
{
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  if (is.null(df)) 
  fac <- qnorm(a)
  else fac <- qt(a, df)
  fac <- c(0, fac)
  fff = format(100 * a, 3)
  pct = paste(fff, "%")
  pct = c("value", pct)
  ci <- array(NA, dim = c(length(parm), 3L), dimnames = list(parm, 
  pct))
  ci[] <- cf + ses %o% fac
  ci
}

# .smartsub <- function (pat, repl, x) 
# {
#   args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
#   names <- if (is.null(names(args))) 
#   character(length(args))
#   else names(args)
#   dovec <- names %in% vectorize.args
#   do.call("mapply", c(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
#   SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
# }


.nlConfint <- function (obj = NULL, texts, level = 0.95, coeff = NULL, Vcov = NULL, 
          df2 = NULL, x = NULL, silent=FALSE) 
{
  
  # obj = NULL; texts = ""; level = 0.95; coeff = NULL; Vcov = NULL; df2 = NULL; x = NULL; silent=FALSE
  if (!is.null(obj)) {
    co = try(coef(obj), silent = TRUE)
    cond = attr(co, "condition")
    if (is.null(coeff) && (is.null(cond))) 
      coeff = co
    vc = try(vcov(obj), silent = TRUE)
    cond2 = attr(vc, "condition")
    if (is.null(Vcov) && (is.null(cond2))) 
      Vcov = vc
  }
  if (is.null(coeff)) {
    if (is.null(obj)) 
      mess = "Both  'obj' and 'coeff' are missing"
    else {
      clm = class(obj)
      part1 = "There are no coef() methods for model objects of class \""
      mess = paste0(part1, clm, "\".\nInput the 'coeff' parameter.")
    }
    stop(mess)
  }
  if (is.null(Vcov)) {
    if (is.null(obj)) 
      mess = "Both  'obj' and 'Vcov' are missing"
    else {
      clm = class(obj)
      part1 = "There are no vcov() methods for model objects of class \""
      mess = paste0(part1, clm, "\".\nInput the 'Vcov' parameter.")
    }
    stop(mess)
  }
  if (length(texts) > 1) kkk = texts[1] else kkk = strsplit(texts[1], ";")[[1]]
  kkkfl = as.formula(paste("~", kkk[1]))
  vvss = setdiff(all.vars(kkkfl), "x")
  # texts = getFromNamespace(".smartsub", ns="nlWaldTest")(vvss, "b", texts)
  if (length(texts) > 1) ltext0 = texts else ltext0 = strsplit(texts, ";")[[1]]
  texts1 = gsub("[", "", texts, fixed = TRUE)
  texts1 = gsub("]", "", texts1, fixed = TRUE)
  if (length(texts1) > 1) ltext = texts1 else ltext = strsplit(texts1, ";")[[1]]
  r = length(ltext)
  n = length(coeff)
  namess = paste0("b", 1:n)
  for (j in 1L:n) assign(namess[j], coeff[j])
  if (!is.null(x)) {
    nx = length(x)
    namesx = paste0("x", 1:nx)
    for (j in 1L:nx) assign(namesx[j], x[j])
  }
  grad = c()
  hess = c()
  for (i in 1L:r) {
    fli <- as.formula(paste("~", ltext[i]))
    z = try(deriv(as.formula(fli), namess), silent = T)
    if (inherits(z, "try-error")) {
      tei = as.character(i)
      tri2 = ", numerical derivatives were used in delta-method"
      wate = paste0("Note: For function ", i, tri2)
      if (!silent) message(wate)
      ez <- numericDeriv(expr=quote(eval(parse(text = ltext[i]))), 
                        theta=namess)
    }
    else ez = eval(z)
    hessj = attr(ez, "gradient")
    grad = rbind(grad, ez[1])
    hess = rbind(hess, hessj)
  }
  Rb = grad
  ddd = hess %*% Vcov %*% t(hess)
  matr = chol2inv(chol(ddd))
  ses = sqrt(diag(ddd))
  trydf = identical(df2, TRUE)
  if (trydf) {
    isdf = try(df.residual(obj), silent = TRUE)
    df2 = isdf
    if (is.null(df2)) {
      wn = "Note: Failed to extract df for denominator; z-intervals applied"
      if (!silent) message(wn)
    }
  }
  getFromNamespace(".getint", ns="HelpersMG")(as.numeric(Rb), ltext0, ses, level, df = df2)
}
