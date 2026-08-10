# Normalizes RDS-loaded objects that may be either a bare Seurat object or
# a list of them, depending on which stage produced the file.
#
# Seurat objects are S4, so naive checks like `class(x) == "list"` or
# `is.list(x)` are not a reliable way to tell "bare object" from "list of
# objects" -- they happen to work for a plain Seurat object (class() and
# is.list() both correctly say "not a list"), but that's incidental, not
# guaranteed, and gives no positive signal about what we actually expect.
# Testing directly for the Seurat class is the only check that's actually
# sensical here.

# Ensure x is a list of Seurat object(s), wrapping a bare object if needed.
as_seurat_list <- function(x) {
  if (inherits(x, "Seurat")) {
    return(list(x))
  }
  x
}

# Ensure x is a single Seurat object, unwrapping a one-(or-more)-element
# list if needed (takes the first element).
first_seurat <- function(x) {
  if (inherits(x, "Seurat")) {
    return(x)
  }
  x[[1]]
}
