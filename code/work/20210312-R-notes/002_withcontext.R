
x <- search()
detach(name = "package:IRanges", character.only = TRUE)
detach(name = "package:AnnotationDbi", character.only = TRUE)
detach("package:org.Mm.eg.db", unload = FALSE)

# install.packages("crunch")
requireNamespace("crunch")


base::detach("package:dplyr")
select

detach_packages <- function(F) {
  # install.packages("crunch")
  requireNamespace("crunch")
  
  base::with(
    crunch::ContextManager(
      enter = function() {
        base::search()
      },
      exit = function() {
        for (pkg in base::rev(base::search())) {
          if (!(pkg %in% enter)) {
            base::detach(pkg, character.only = TRUE)
          }
        }
      },
      # this overwrites the variable `enter`, if existent
      as = "enter"
    ),
    {
      F()
    }
  )
}
