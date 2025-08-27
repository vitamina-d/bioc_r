#* echo back the input
#* @param msg The message to echo
#* @get /echo
#* @serializer unboxedJSON 
function(msg = "msg") {
  list(
    msg = paste0("The message is: '", msg, "'")
  )
}