# Shiny App for Molecular Weight Calculation and Lipinski's Rule Check
# Author: Karol Topka Klonczynski
# Date: 2023-05-04

# Required libraries
library(shiny)
library(stringr)
library(httr)

# Atomic mass dictionary for elements
masses <- list(
  "Ru" = "101.07(2)", "Re" = "186.207(1)", "Rf" = "[267]", "Rg" = "[280]",
  "Ra" = "[226]", "Rb" = "85.4678(3)", "Rn" = "[222]", "Rh" = "102.90550(2)",
  "Be" = "9.012182(3)", "Ba" = "137.327(7)", "Bh" = "[272]", "Bi" = "208.98040(1)",
  "Bk" = "[247]", "Br" = "79.904(1)", "Uuh" = "[293]", "H" = "1.00794(4)",
  "P" = "30.973762(2)", "Os" = "190.23(3)", "Es" = "[252]", "Hg" = "200.59(2)",
  "Ge" = "72.64(1)", "Gd" = "157.25(3)", "Ga" = "69.723(1)", "Pr" = "140.90765(2)",
  "Pt" = "195.084(9)", "Pu" = "[244]", "C" = "12.0107(8)", "Pb" = "207.2(1)",
  "Pa" = "231.03588(2)", "Pd" = "106.42(1)", "Cd" = "112.411(8)", "Po" = "[209]",
  "Pm" = "[145]", "Hs" = "[270]", "Uuq" = "[289]", "Uup" = "[288]", "Uus" = "",
  "Uuo" = "[294]", "Ho" = "164.93032(2)", "Hf" = "178.49(2)", "K" = "39.0983(1)",
  "He" = "4.002602(2)", "Md" = "[258]", "Mg" = "24.3050(6)", "Mo" = "95.96(2)",
  "Mn" = "54.938045(5)", "O" = "15.9994(3)", "Mt" = "[276]", "S" = "32.065(5)",
  "W" = "183.84(1)", "Zn" = "65.38(2)", "Eu" = "151.964(1)", "Zr" = "91.224(2)",
  "Er" = "167.259(3)", "Ni" = "58.6934(4)", "No" = "[259]", "Na" = "22.98976928(2)",
  "Nb" = "92.90638(2)", "Nd" = "144.242(3)", "Ne" = "20.1797(6)", "Np" = "[237]",
  "Fr" = "[223]", "Fe" = "55.845(2)", "Fm" = "[257]", "B" = "10.811(7)",
  "F" = "18.9984032(5)", "Sr" = "87.62(1)", "N" = "14.0067(2)", "Kr" = "83.798(2)",
  "Si" = "28.0855(3)", "Sn" = "118.710(7)", "Sm" = "150.36(2)", "V" = "50.9415(1)",
  "Sc" = "44.955912(6)", "Sb" = "121.760(1)", "Sg" = "[271]", "Se" = "78.96(3)",
  "Co" = "58.933195(5)", "Cn" = "[285]", "Cm" = "[247]", "Cl" = "35.453(2)",
  "Ca" = "40.078(4)", "Cf" = "[251]", "Ce" = "140.116(1)", "Xe" = "131.293(6)",
  "Lu" = "174.9668(1)", "Cs" = "132.9054519(2)", "Cr" = "51.9961(6)", "Cu" = "63.546(3)",
  "La" = "138.90547(7)", "Li" = "6.941(2)", "Tl" = "204.3833(2)", "Tm" = "168.93421(2)",
  "Lr" = "[262]", "Th" = "232.03806(2)", "Ti" = "47.867(1)", "Te" = "127.60(3)",
  "Tb" = "158.92535(2)", "Tc" = "[98]", "Ta" = "180.94788(2)", "Yb" = "173.054(5)",
  "Db" = "[268]", "Dy" = "162.500(1)", "Ds" = "[281]", "I" = "126.90447(3)",
  "U" = "238.02891(3)", "Y" = "88.90585(2)", "Ac" = "[227]", "Ag" = "107.8682(2)",
  "Uut" = "[284]", "Ir" = "192.217(3)", "Am" = "[243]", "Al" = "26.9815386(8)",
  "As" = "74.92160(2)", "Ar" = "39.948(1)", "Au" = "196.966569(4)", "At" = "[210]",
  "In" = "114.818(3)"
)

# Function to tokenize the molecule formula
tokenize_molecule <- function(formula) {
  # Tokenize molecule formula into elements and numbers
  regmatches(formula, gregexpr("([A-Z][a-z]?|\\d+|\\(|\\)|\\[|\\])", formula, perl=TRUE))[[1]]
}

# Function to parse the formula and calculate element counts
parse_formula <- function(tokens) {
  stack <- list()
  elements <- list()
  
  i <- 1
  while (i <= length(tokens)) {
    token <- tokens[i]
    
    if (str_detect(token, "^[A-Z]")) {
      if (i + 1 <= length(tokens) && str_detect(tokens[i + 1], "^\\d+$")) {
        count <- as.integer(tokens[i + 1])
        i <- i + 1
      } else {
        count <- 1
      }
      
      # Increment the count for the element
      if (token %in% names(elements)) {
        elements[[token]] <- elements[[token]] + count
      } else {
        elements[[token]] <- count
      }
    } else if (token == "(" || token == "[") {
      stack <- c(stack, list(elements))
      elements <- list()
    } else if (token == ")" || token == "]") {
      subelements <- elements
      elements <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      
      if (i + 1 <= length(tokens) && str_detect(tokens[i + 1], "^\\d+$")) {
        repeat_count <- as.integer(tokens[i + 1])
        i <- i + 1
      } else {
        repeat_count <- 1
      }
      
      # Multiply the counts of subelements by repeat_count and add to elements
      for (element in names(subelements)) {
        if (element %in% names(elements)) {
          elements[[element]] <- elements[[element]] + subelements[[element]] * repeat_count
        } else {
          elements[[element]] <- subelements[[element]] * repeat_count
        }
      }
    }
    
    i <- i + 1
  }
  
  return(elements)
}

# Function to calculate the molecular weight and handle uncertain elements
calculate_mass <- function(formula, masses) {
  tokens <- tokenize_molecule(formula)
  elements <- parse_formula(tokens)
  
  total_mass <- 0.0
  uncertain_elements <- c()
  
  result <- data.frame(Element = character(), Count = integer(), Contribution = character(), stringsAsFactors = FALSE)
  
  for (element in names(elements)) {
    count <- elements[[element]]
    if (element %in% names(masses)) {
      mass_str <- masses[[element]]
      if (startsWith(mass_str, "[")) {
        uncertain_elements <- c(uncertain_elements, element)
        result <- rbind(result, data.frame(Element = element, Count = count, Contribution = "[Uncertain]"))
      } else {
        mass <- as.numeric(str_replace(mass_str, "\\(.*\\)", ""))
        contribution <- mass * count
        total_mass <- total_mass + contribution
        result <- rbind(result, data.frame(Element = element, Count = count, Contribution = sprintf("%.4f u", contribution)))
      }
    } else {
      warning(sprintf("Element '%s' not found in mass dictionary.", element))
    }
  }
  
  cat("Detailed Calculation:\n")
  print(result, row.names = FALSE)
  cat(sprintf("\nTotal mass: %.6f u\n", total_mass))
  
  invisible(list(total_mass = total_mass, uncertain_elements = uncertain_elements))
}

check_lipinski_rule <- function(formula) {
  tokens <- tokenize_molecule(formula)
  elements <- parse_formula(tokens)
  
  molecular_weight <- 0
  n_atoms <- 0
  o_atoms <- 0
  
  for (element in names(elements)) {
    count <- elements[[element]]
    if (element %in% names(masses)) {
      mass <- as.numeric(str_replace(masses[[element]], "\\(.*\\)", ""))
      molecular_weight <- molecular_weight + mass * count
    }
    if (element == "N") {
      n_atoms <- count
    }
    if (element == "O") {
      o_atoms <- count
    }
  }
  
  mw_rule <- molecular_weight <= 500
  no_rule <- (n_atoms + o_atoms) <= 10
  
  passes_lipinski <- mw_rule && no_rule
  
  return(list(
    passes_lipinski = passes_lipinski,
    mw_rule = mw_rule,
    no_rule = no_rule
  ))
}
# Shiny UI
ui <- fluidPage(
  titlePanel("Molecular Weight Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("formula", "Chemical Formula:", value = ""),
      actionButton("calculate", "Calculate")
    ),
    
    mainPanel(
      verbatimTextOutput("result"),
      h4("Lipinski's Rule (without logP and NH and OH groups rule):"),
      textOutput("lipinski_result"),
      textOutput("mw_result"),
      textOutput("no_result"),
      textOutput("nhoh_result"),
      h4("Recently calculated compounds:"),
      uiOutput("recent_formulas")
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  recent_formulas <- reactiveVal(character())
  
  observeEvent(input$calculate, {
    formula <- input$formula
    result <- capture.output(calculate_mass(formula, masses))
    output$result <- renderText({
      paste(result, collapse = "\n")
    })
    
    # Check Lipinski's rule (without logP)
    lipinski_results <- check_lipinski_rule(formula)
    
    output$lipinski_result <- renderText({
      if (!is.null(lipinski_results$passes_lipinski)) {
        if (lipinski_results$passes_lipinski) {
          "The compound passes Lipinski's Rule (without logP and NH and OH groups rule)."
        } else {
          "The compound does not pass Lipinski's Rule (without logP and NH and OH groups rule)."
        }
      } else {
        "Unable to assess Lipinski's Rule."
      }
    })
    
    output$mw_result <- renderText({
      if (!is.null(lipinski_results$mw_rule)) {
        if (lipinski_results$mw_rule) {
          "Molecular Weight (<=500): PASS"
        } else {
          "Molecular Weight (<=500): FAIL"
        }
      } else {
        "Unable to assess Molecular Weight rule."
      }
    })
    
    output$no_result <- renderText({
      if (!is.null(lipinski_results$no_rule)) {
        if (lipinski_results$no_rule) {
          "Number of N and O atoms (<=10): PASS"
        } else {
          "Number of N and O atoms (<=10): FAIL"
        }
      } else {
        "Unable to assess N and O atoms rule."
      }
    })
    
    recent_formulas(c(formula, recent_formulas()))
  })
  updateFormulaInput <- function(formula) {
    updateTextInput(session, "formula", value = formula)
  }
  
  output$recent_formulas <- renderUI({
    formulas <- recent_formulas()
    if (length(formulas) > 0) {
      tags$ul(
        lapply(formulas, function(formula) {
          tags$li(
            tags$a(href = "#", onclick = sprintf("Shiny.onInputChange('update_formula', '%s')", formula), formula)
          )
        })
      )
    } else {
      tags$p("No recently calculated compounds.")
    }
  })
  observeEvent(input$update_formula, {
    updateFormulaInput(input$update_formula)
  })
}


# Run the Shiny app
shinyApp(ui, server)