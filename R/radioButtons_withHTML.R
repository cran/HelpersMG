
.radioButtons_withHTML <- function (inputId, label, choices, selected = NULL, 
								    inline = FALSE, width = NULL) 
{

  # http://stackoverflow.com/users/4474157/nice
  
  tags <- NULL
  HTML <- NULL
  div <- NULL
  validateCssUnit <- NULL
  generateOptions_withHTML <- function (inputId, choices, selected, inline, type = "checkbox") 
  {
    options <- mapply(choices, names(choices), FUN = function(value, 
                                                              name) {
      inputTag <- tags$input(type = type, name = inputId, value = value)
      if (value %in% selected) 
        inputTag$attribs$checked <- "checked"
      if (inline) {
        tags$label(class = paste0(type, "-inline"), inputTag, 
                   tags$span(HTML(name)))
      }
      else {
        tags$div(class = type, tags$label(inputTag, tags$span(HTML(name))))
      }
    }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    div(class = "shiny-options-group", options)
  }
  choices <- getFromNamespace("choicesWithNames", ns="shiny")(choices)
    # shiny:::choicesWithNames
  selected <- if (is.null(selected)) 
    choices[[1]]
  else {
    getFromNamespace("validateSelected", ns="shiny")(selected, choices, inputId)
    # shiny:::validateSelected(selected, choices, inputId)
  }
  if (length(selected) > 1) 
    stop("The 'selected' argument must be of length 1")
  options <- generateOptions_withHTML(inputId, choices, selected, inline, 
                                      type = "radio")
  divClass <- "form-group shiny-input-radiogroup shiny-input-container"
  if (inline) 
    divClass <- paste(divClass, "shiny-input-container-inline")
  tags$div(id = inputId, style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"), class = divClass, 
    getFromNamespace("controlLabel", ns="shiny")(inputId, label), options)
}
