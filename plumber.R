library(plumber)

# Create an initial router
route <- pr() |>
  pr_get("/foo", function() "foo")

#* @plumber
function(pr) {
  pr |>
    pr_mount("/bar", route)
}
