#* echo back the input
#* @param msg The message to echo
#* @get /
function(msg = "") {
  list(
    msg = paste0("The message is: '", msg, "'")
  )
}