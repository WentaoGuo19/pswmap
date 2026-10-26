# Argument checks (internal). Each check stops with a message naming the argument, or returns its value.

stop_arg = function(...) stop(..., call. = FALSE)
warn_arg = function(...) warning(..., call. = FALSE)

is_number = function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x)

# One of a set of character choices; the full vector of choices (the default) gives the first one. Partial matching
# as in match.arg().
check_choice = function(x, choices, name) {
  if (identical(x, choices)) return(choices[1])
  i = if (is.character(x) && length(x) == 1 && !is.na(x)) pmatch(x, choices) else NA
  if (is.na(i)) stop_arg("`", name, "` must be one of ", paste0("\"", choices, "\"", collapse = ", "), ".")
  choices[i]
}

check_string = function(x, name) {
  if (!(is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)))
    stop_arg("`", name, "` must be a single column name (a character string).")
  x
}

check_flag = function(x, name) {
  if (!(is.logical(x) && length(x) == 1 && !is.na(x))) stop_arg("`", name, "` must be TRUE or FALSE.")
  x
}

check_number = function(x, name) {
  if (!is_number(x)) stop_arg("`", name, "` must be a single finite number.")
  x
}

check_positive = function(x, name, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(x)
  if (!(is_number(x) && x > 0))
    stop_arg("`", name, "` must be ", if (null_ok) "NULL or ", "a single positive number.")
  x
}

# A number in the open interval (lower, upper), or the half-open [lower, upper) with closed_lower = TRUE.
check_between = function(x, name, lower = 0, upper = 1, closed_lower = FALSE, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(x)
  ok = is_number(x) && (if (closed_lower) x >= lower else x > lower) && x < upper
  if (!ok) stop_arg("`", name, "` must be ", if (null_ok) "NULL or ", "a single number in ",
                    if (closed_lower) "[" else "(", lower, ", ", upper, ").")
  x
}

check_whole = function(x, name, min = 1) {
  if (!(is_number(x) && x == round(x) && x >= min)) stop_arg("`", name, "` must be a whole number of at least ", min, ".")
  as.integer(x)
}

# A list of named arguments (empty by default).
check_arg_list = function(x, name) {
  if (is.null(x)) return(list())
  if (!is.list(x) || is.data.frame(x)) stop_arg("`", name, "` must be a list, e.g. list(name = value).")
  if (length(x) && (is.null(names(x)) || any(names(x) == "") || anyDuplicated(names(x))))
    stop_arg("every element of `", name, "` must have a distinct name, e.g. list(name = value).")
  x
}

# An RBesT normal mixture, e.g. from RBesT::mixnorm().
check_normmix = function(x, name, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(x)
  if (!inherits(x, "normMix"))
    stop_arg("`", name, "` must be ", if (null_ok) "NULL or ", "an RBesT normal mixture, e.g. RBesT::mixnorm(c(1, 0, 10), sigma = 1).")
  x
}

# Outcomes of one trial arm: at least two finite values that are not all equal.
check_arm = function(y, name) {
  if (!is.numeric(y) || length(y) < 2 || anyNA(y) || any(!is.finite(y)))
    stop_arg("`", name, "` must be a numeric vector of at least two finite outcomes.")
  if (stats::var(y) == 0) stop_arg("the outcomes in `", name, "` are all equal.")
  y
}

# Arguments the caller typed (names of match.call(), without the function name).
given_args = function(call) setdiff(names(as.list(call))[-1], "")

# Warn that arguments the caller typed have no effect.
warn_ignored = function(args, given, reason) {
  used = intersect(args, given)
  if (length(used))
    warn_arg(paste0("`", used, "`", collapse = ", "), if (length(used) == 1) " is" else " are", " ignored ", reason, ".")
}
