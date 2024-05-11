# Shiny App for Molecular Weight Calculation and Lipinski's Rule Check
# Author: Karol Topka Klonczynski
# Date: 2023-05-04

# Required libraries
library(shiny)
library(stringr)
library(httr)
library(jsonlite)
library(rcdk)

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

# Initialize an empty list to store element multiplicities
multiples <- list()

# Function to add commas between element symbols
add_commas <- function(orig) {
  with_commas <- str_replace_all(orig, "(?<=\\p{Lu})(?=\\p{Lu})", ",")
  paste0(with_commas, ",")
}

# Function to find the number of occurrences of an element
find_num <- function(i, val) {
  if (str_detect(str_sub(val, i, i), "[A-Za-z,]")) {
    return(str_sub(val, i, i))
  }
  count <- i
  while (count < nchar(val) && !str_detect(str_sub(val, count, count), "[A-Za-z,]")) {
    count <- count + 1
  }
  if (count >= nchar(val)) {
    return(str_sub(val, i))
  }
  str_sub(val, i, count - 1)
}

# Function to distribute the number of occurrences of an element
distribute <- function(val, num) {
  build <- ""
  sub <- 1
  for (i in seq_len(nchar(val))) {
    if (!str_detect(str_sub(val, i, i), "[A-Za-z]")) {
      sub <- as.numeric(find_num(i, val))
      sub <- sub * num
      build <- paste0(build, sub)
      i <- i + nchar(find_num(i, val)) - 1
    } else if (i + 1 >= nchar(val) || str_detect(str_sub(val, i + 1, i + 1), "\\p{Lu}")) {
      sub <- 1 * num
      build <- paste0(build, str_sub(val, i, i), sub)
      i <- i + nchar(find_num(i, val))
      sub <- 1
    } else {
      build <- paste0(build, str_sub(val, i, i))
    }
  }
  build
}

# Function to distribute the number of occurrences in parentheses
distrib_parenth <- function(val) {
  op <- str_locate(val, "\\(")[1]
  num <- 1
  while (!is.na(op)) {
    close <- str_locate(val, "\\)", op)[1]
    between <- str_sub(val, op + 1, close - 1)
    if (close + 1 < nchar(val) && !str_detect(str_sub(val, close + 1, close + 1), "[A-Za-z]") && str_sub(val, close + 1, close + 1) != "") {
      num <- as.numeric(str_sub(val, close + 1, close + 1))
      val <- paste0(str_sub(val, 1, op - 1), distribute(between, num), str_sub(val, close + 2))
    } else {
      val <- paste0(str_sub(val, 1, op - 1), distribute(between, num), str_sub(val, close + 1))
    }
    op <- str_locate(val, "\\(", op + 1)[1]
    num <- 1
  }
  val
}

# Function to add markers to the formula
add_markers <- function(val) {
  if (str_detect(val, "\\[") && str_detect(val, "\\]")) {
    val <- str_replace_all(val, "\\[", "(")
    val <- str_replace_all(val, "\\]", ")")
    val <- distrib_parenth(val)
  }
  if (!str_detect(val, ",")) {
    val <- add_commas(val)
  }
  val
}

# Function to add element multiplicities to the 'multiples' list
add <- function(mult, symb) {
  com <- str_locate(symb, ",")[1]
  if (!is.na(com)) {
    symb <- str_sub(symb, 1, com - 1)
  }
  if (!(symb %in% names(multiples))) {
    multiples[[symb]] <<- mult
  } else {
    multiples[[symb]] <<- multiples[[symb]] + mult
  }
}

# Function to strip coefficients from the formula
strip_coeff <- function(val) {
  build <- ""
  symb <- ""
  coeff <- FALSE
  i <- 1
  while (i <= nchar(val)) {
    if (!str_detect(str_sub(val, i, i), "[A-Za-z]") && str_sub(val, i, i) != ",") {
      num <- as.numeric(find_num(i, val))
      add(num, symb)
      coeff <- TRUE
    }
    if (str_sub(val, i, i) == ",") {
      if (str_detect(str_sub(val, i, i), "[A-Za-z]") || str_sub(val, i, i) == ",") {
        symb <- paste0(symb, str_sub(val, i, i))
      }
      if (!coeff) {
        add(1, symb)
      }
      build <- paste0(build, symb)
      symb <- ""
      coeff <- FALSE
    } else if (str_detect(str_sub(val, i, i), "[A-Za-z]")) {
      symb <- paste0(symb, str_sub(val, i, i))
    }
    i <- i + nchar(find_num(i, val))
  }
  if (nchar(build) == 0) {
    return("")
  }
  build
}

# Function to calculate the mass of a single element
calc_single_element <- function(val) {
  count <- 0
  num <- as.numeric(str_replace(masses[[val]], " g/mol", ""))
  if (val %in% names(multiples)) {
    count <- count + num * multiples[[val]]
  } else {
    count <- count + num
  }
  count
}

# Function to calculate the total molecular weight
calculate <- function(comp) {
  build <- ""
  count <- 0
  finished <- c()
  for (i in seq_len(nchar(comp))) {
    if (str_sub(comp, i, i) == ",") {
      if (build %in% names(masses) && !(build %in% finished)) {
        count <- count + calc_single_element(build)
        finished <- c(finished, build)
      } else if (!(build %in% names(masses))) {
        return(build)
      }
      build <- ""
    } else if (i == nchar(comp)) {
      build <- paste0(build, str_sub(comp, i, i))
      if (build %in% names(masses) && !(build %in% finished)) {
        count <- count + calc_single_element(build)
        finished <- c(finished, build)
      } else if (!(build %in% names(masses))) {
        return(build)
      }
    } else {
      build <- paste0(build, str_sub(comp, i, i))
    }
  }
  count
}

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

# Function to check Lipinski's Rule (without logP)
check_lipinski_rule <- function(formula) {
  n_atoms <- str_count(formula, "N")
  o_atoms <- str_count(formula, "O")
  nh_groups <- n_atoms
  oh_groups <- o_atoms
  molecular_weight <- calculate_mass(formula, masses)$total_mass
  mw_rule <- molecular_weight <= 500
  no_rule <- (n_atoms + o_atoms) < 10
  nhoh_rule <- (nh_groups + oh_groups) <= 5
  passes_lipinski <- mw_rule && no_rule && nhoh_rule
  return(list(
    passes_lipinski = passes_lipinski,
    mw_rule = mw_rule,
    no_rule = no_rule,
    nhoh_rule = nhoh_rule
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
      h4("Lipinski's Rule (without logP):"),
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
      if (lipinski_results$passes_lipinski) {
        "The compound passes Lipinski's Rule (without logP)."
      } else {
        "The compound does not pass Lipinski's Rule (without logP)."
      }
    })
    
    output$mw_result <- renderText({
      if (lipinski_results$mw_rule) {
        "Molecular Weight (<=500): PASS"
      } else {
        "Molecular Weight (<=500): FAIL"
      }
    })
    
    output$no_result <- renderText({
      if (lipinski_results$no_rule) {
        "Number of N and O atoms (<10): PASS"
      } else {
        "Number of N and O atoms (<10): FAIL"
      }
    })
    
    output$nhoh_result <- renderText({
      if (lipinski_results$nhoh_rule) {
        "Number of NH and OH groups (<=5): PASS"
      } else {
        "Number of NH and OH groups (<=5): FAIL"
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