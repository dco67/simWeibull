library(shiny)
library(shinyWidgets)
library(ggplot2)
library(DT)
library(fitdistrplus)
library(kableExtra)
library(dplyr)
library(ForestFit)
library(EnvStats)

weibull_characteristics <- function(shape, scale) {
  mean_weibull <- scale * gamma(1 + 1/shape)
  variance_weibull <- (scale^2) * (gamma(1 + 2/shape) - (gamma(1 + 1/shape))^2)
  sd_weibull <- sqrt(variance_weibull)
  skewness_weibull <- (gamma(1 + 3/shape) - 3*mean_weibull*gamma(1 + 2/shape) + 2*mean_weibull^3) / (sd_weibull^3)
  kurtosis_weibull <- (-6*(gamma(1 + 1/shape))^2 - 3*(gamma(1 + 2/shape))^2 + gamma(1 + 4/shape)) / (variance_weibull^2)
  median_weibull <- scale * (log(2))^(1/shape)
  mode_weibull <- if (shape > 1) scale * ((shape - 1)/shape)^(1/shape) else 0
  entropy_weibull <- (1 - 1/shape) + log(scale) + gamma(1 + 1/shape) * (1 - 1/shape)
  moment1 <- mean_weibull
  moment2 <- scale^2 * gamma(1 + 2/shape)
  moment3 <- scale^3 * gamma(1 + 3/shape)
  moment4 <- scale^4 * gamma(1 + 4/shape)
  result <- data.frame(
    Characteristic = c("Mean", "Variance", "Standard Deviation", "Skewness", "Kurtosis", "Median", "Mode", "Entropy", "1st Moment", "2nd Moment", "3rd Moment", "4th Moment"),
    Value = c(mean_weibull, variance_weibull, sd_weibull, skewness_weibull, kurtosis_weibull, median_weibull, mode_weibull, entropy_weibull, moment1, moment2, moment3, moment4)
  )
  return(result)
}

dweibull3 <- function(x, mu, beta, lambda) {
  ifelse(x < mu, 0, (beta / lambda) * ((x - mu) / lambda)^(beta - 1) * exp(-((x - mu) / lambda)^beta))
}

pweibull3 <- function(x, mu, beta, lambda) {
  ifelse(x < mu, 0, 1 - exp(-((x - mu) / lambda)^beta))
}

qweibull3 <- function(p, mu, beta, lambda) {
  mu + lambda * (-log(1 - p))^(1 / beta)
}

rweibull3 <- function(n, mu, beta, lambda) {
  u <- runif(n)  # Generate n uniform random numbers between 0 and 1
  mu + lambda * (-log(1 - u))^(1 / beta)  # Apply the inverse CDF transformation
}


# Define the table data with equations
equation_table <- data.frame(
  Characteristic = c("Moyenne", "Variance", "Écart-type", "Asymétrie", "Voussure", "Médiane", "Mode", "Entropie", "1er Moment", "2ème Moment", "3ème Moment", "4ème Moment"),
  Equation = c(
    "$\\mu = \\beta \\cdot \\Gamma\\left(1 + \\frac{1}{\\alpha}\\right)$",
    "$\\sigma^2 = \\beta^2 \\cdot \\left[ \\Gamma\\left(1 + \\frac{2}{\\alpha}\\right) - \\left(\\Gamma\\left(1 + \\frac{1}{\\alpha}\\right)\\right)^2 \\right]$",
    "$\\sigma = \\sqrt{\\sigma^2}$",
    "$\\gamma_1 = \\frac{\\Gamma\\left(1 + \\frac{3}{\\alpha}\\right) - 3\\mu\\Gamma\\left(1 + \\frac{2}{\\alpha}\\right) + 2\\mu^3}{\\sigma^3}$",
    "$\\gamma_2 = \\frac{-6\\left(\\Gamma\\left(1 + \\frac{1}{\\alpha}\\right)\\right)^2 - 3\\left(\\Gamma\\left(1 + \\frac{2}{\\alpha}\\right)\\right)^2 + \\Gamma\\left(1 + \\frac{4}{\\alpha}\\right)}{\\sigma^4}$",
    "$\\tilde{x} = \\beta \\cdot \\left(\\log(2)\\right)^{\\frac{1}{\\alpha}}$",
    "$M = \\begin{cases} \\beta \\cdot \\left(\\frac{\\alpha - 1}{\\alpha}\\right)^{\\frac{1}{\\alpha}} & \\text{si } \\alpha > 1 \\\\ 0 & \\text{sinon} \\end{cases}$",
    "$H = 1 - \\frac{1}{\\alpha} + \\log(\\beta) + \\Gamma\\left(1 + \\frac{1}{\\alpha}\\right)\\left(1 - \\frac{1}{\\alpha}\\right)$",
    "$m_1 = \\mu$",
    "$m_2 = \\beta^2 \\cdot \\Gamma\\left(1 + \\frac{2}{\\alpha}\\right)$",
    "$m_3 = \\beta^3 \\cdot \\Gamma\\left(1 + \\frac{3}{\\alpha}\\right)$",
    "$m_4 = \\beta^4 \\cdot \\Gamma\\left(1 + \\frac{4}{\\alpha}\\right)$"
  )
)

equation_table_html <- equation_table %>%
  kbl(format = "html", escape = FALSE, col.names = c("Caractéristique", "Équation")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width=FALSE) %>%
  column_spec(1, bold=TRUE, border_left=TRUE) %>%
  column_spec(2, bold=FALSE, border_right=TRUE) %>%
  row_spec(12, extra_css = "border-bottom: 1px solid;")


# References:
# https://statisticsbyjim.com/probability/weibull-distribution/
# https://stats.libretexts.org/Bookshelves/Probability_Theory/Probability_Mathematical_Statistics_and_Stochastic_Processes_(Siegrist)/05%3A_Special_Distributions/5.38%3A_The_Weibull_Distribution

ui <- fluidPage(
  titlePanel(h3("La distribution de Weibull")),
  withMathJax(),
  tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(
          6,
          switchInput(
            inputId = "input_type",
            label = "Input",
            onLabel = "Slider",
            offLabel = "Numeric",
            value = TRUE,
            size = "mini"
          )
        ),
        column(4,
               offset = 2,
               actionBttn("reset_button",
                          "Réinit",
                          icon = icon("refresh"),
                          class = "btn-danger",
                          style = "bordered",
                          size = "xs",
                          color = "default"
               )
        )
      ),
      conditionalPanel(
        condition = "input.input_type == true",
        sliderInput("shape",
                    "Forme [$\\alpha$]:",
                    min = 0.1,
                    max = 100,
                    step = 0.1,
                    value = 2
        ),
        sliderInput("scale",
                    "Échelle [$\\beta$]:",
                    min = 0.1,
                    max = 100,
                    step = 0.1,
                    value = 1
        )
      ),
      conditionalPanel(
        condition = "input.input_type == false",
        numericInput("shape_num",
                     "Forme [$\\alpha$]:",
                     min = 0.1,
                     step = 0.1,
                     value = 2,
                     width = 150
        ),
        numericInput("scale_num",
                     "Échelle [$\\beta$]:",
                     min = 0.1,
                     step = 0.1,
                     value = 1,
                     width = 150
        )
      ),
      switchInput(
        inputId = "input_mode",
        label = "Data",
        onLabel = "Générées",
        offLabel = "Téléchargées",
        value = TRUE,
        size = "mini"
      ),
      conditionalPanel(
        condition = "input.input_mode == true",
        numericInput("num_observations",
                     "# d'observations:",
                     value = 1000,
                     min = 1,
                     max = 10000,
                     width = 150
        ),
        fluidRow(
          column(
            6,
            actionBttn(
              inputId = "generate_button",
              label = "Générer",
              style = "unite",
              color = "primary",
              icon = icon("sliders"),
              size = "sm"
            )
          ),
          column(
            6,
            downloadBttn(
              outputId = "download_data",
              label = "Télécharger",
              style = "unite",
              color = "default",
              icon = icon("download"),
              size = "sm"
            )
          )
        )
      ),
      conditionalPanel(
        condition = "input.input_mode == false",
        fileInput("file1", "File:",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                  )
        ),
        DTOutput("fittedParams")
      ),
      p(),
      wellPanel(style = "background: lightblue",
                fluidRow(
                  column(4,
                         a(h4("Par Daniel Coulombe, Ph.D.")),
                         p("2024")
                  ),
                  column(4,
                         tags$a(
                           href="https://isteah.org", 
                           tags$img(src="ISTEAH_LOGO.png", 
                                    title="ISTEAH", 
                                    width="160",
                                    height="140")
                           )
                         )
                  )
                )
      ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          strong(h4("Introduction")),
          helpText(
            "La distribution de Weibull est une distribution de probabilité continue largement utilisée en ingénierie, entre autre. Par exemple, on utilisera cette distribution ",
            "pour analyser la fiabilité d'un produit ou d'un instrument, ou pour examiner la durée de vie de ces produits et leurs temps de défaillance. Elle trouve des applications dans plusieurs ",
            "autres champs d'intérêt scientifique. ", p(),
            strong("Définition"), p(),
            "Elle est définie principalement par deux paramètres : la forme (\\alpha) et l'échelle (\\beta).", br(),
            p("La fonction de densité de probabilité [PDF] de la distribution de Weibull est donnée par :"),
            tags$p(HTML("$$f(x; \\alpha, \\beta) = \\frac{\\alpha}{\\beta} \\left( \\frac{x}{\\beta} \\right)^{\\alpha-1} e^{-\\left( \\frac{x}{\\beta} \\right)^\\alpha}$$")),
            "La fonction de distribution cumulative [CDF], que l'on utilisera pour déterminer la proportion des observations inférieures à une valeur spécifiée, est définie par :", p(),
            tags$p(HTML("$$F(x; \\alpha, \\beta) = 1 - e^{-\\left(\\frac{x}{\\beta}\\right)^\\alpha}$$")), p(),
            "et la fonction quantile, permettant d'obtenir la valeur sous laquelle on trouve une proportion donnée des observations, est: ", p(),
            tags$p(HTML("$$Q(p; \\alpha, \\beta) = \\alpha \\left[-ln(1-p) \\right ]^{\\frac{1}{\\beta}} $$")), p(),
            p("où :"),
              "$\\alpha$ est le paramètre de forme.",  br(),
              "$\\beta$ est le paramètre d'échelle.",  p(),
            
            "Il est souvent utile d'ajouter un troisième paramètre, $\\gamma$, indiquant la location ou le seuil de la distribution, sur l'abcisse. Dans ce cas, on définit: ", p(),
            tags$p(HTML("$$f(x; \\alpha, \\beta, \\gamma) = \\frac{\\alpha}{\\beta} \\left( \\frac{x-\\gamma}{\\beta} \\right)^{\\alpha-1} e^{-\\left( \\frac{x-\\gamma}{\\beta} \\right)^\\alpha}$$")),

            "la fonction cumulative est définie par: ", p(),
            tags$p(HTML("$$F(x; \\alpha, \\beta, \\gamma) = 1 - e^{-\\left(\\frac{x-\\gamma}{\\beta}\\right)^\\alpha}$$")), p(),
            
            "et la fonction quantile est:", p(),
            tags$p(HTML("$$Q(p; \\alpha, \\beta, \\gamma) = \\gamma+\\alpha \\left[-ln(1-p) \\right ]^{\\frac{1}{\\beta}} $$")), p(),
            p(),
            
            strong("Formes générales"), p(),
            "1. Pour $0<\\alpha<1$, $f(x)$ est décroissante et convexe, avec un mode égal à $\\displaystyle{\\lim_{x \\to 0}} f(x)=\\infty$;", p(),
            "2. Pour $\\alpha = 1$, $f(x)$ est décroissante et convexe, avec un mode nul;", p(),
            "3. Pour $\\alpha>1$, $f(x)$ est concave [croissante puis décroissante], avec un mode égal à $\\left(\\frac{\\alpha-1}{\\alpha}\\right)^{\\frac{1}{\\alpha}}$;", p(),
            "4. Pour $1<\\alpha\\le 2$, $f(x)$ est concave puis convexe, avec un point d'inflexion à $x = \\left[\\frac{3 (\\alpha - 1) + \\sqrt{(5 \\alpha - 1)(\\alpha - 1)}}{2 \\alpha}\\right]^{1/\\alpha}$", p(),
            
           "Les graphiques suivants illustrent différentes formes que peut emprunter une distribution de Weibull, en fonction de ses paramètres: "  , p(),
           p(),
           
            fluidRow(
              column(6,
                     plotOutput("shapePlot")),
              column(6,
                     plotOutput("scalePlot"))
            ),
            p(),
           fluidRow(
             column(6, offset = 3,
                    plotOutput("locationPlot")
                    )
             ),
           
           strong("Caractéristiques"), p(),

           p("Les caractéristiques de la distribution de Weibull sont résumées dans le tableau ci-dessous :"),
           HTML(equation_table_html),
           
           p(),
           "où $\\Gamma(x)$ est la fonction gamma. Rappelons que cette fonction est une généralisation de la fonction factorielle: ", p(),
           "$$\\Gamma(x)= (x-1)! $$",
           
           "De manière générale, pour tout nombre réel $x$, on trouve:", p(),
           "$$\\Gamma(x) = \\int_0^\\infty t^{x-1} e^{-t} \\, dt $$",
           p(),
           
           strong("Domaines d'application"), p(),
           
            "La distribution de Weibull peut modéliser différents types de taux de défaillance, y compris des taux de défaillance croissants, constants et décroissants.", p(),
            "En ingénierie de la fiabilité, le facteur de forme mesure la forme de la courbe de taux de défaillance au fil du temps. ",
            "Le taux de défaillance est le taux auquel un dispositif ou un système échoue sur une période donnée. Il est généralement exprimé comme ",
            "le nombre de défaillances par unité de temps, tel que les défaillances par heure ou les défaillances par année.", p(),
            "Si le facteur de forme est inférieur à 1, la courbe de taux de défaillance est dite décroissante ou diminuante. ",
            "Cela signifie que le taux de défaillances diminue au fil du temps. Cela peut se produire si le dispositif ou le système est bien conçu ",
            "et a peu de défaillances initiales ou si le dispositif ou le système est réparé ou remplacé après une défaillance.", p(),  
            "Si le facteur de forme est supérieur à 1, la courbe de taux de défaillance est dite croissante ou croissante. ",
            "Cela signifie que le taux de défaillances augmente au fil du temps. Cela peut se produire si le dispositif ou le système n'est pas bien entretenu ",
            "ou s'il est soumis à des conditions extrêmes qui augmentent la probabilité de défaillance.", p(),
            "Si le facteur de forme est égal à 1, la courbe de taux de défaillance est dite constante. Cela signifie que le taux de défaillances est constant ",
            "au fil du temps. Cela peut se produire si le dispositif ou le système est bien conçu et est soumis à des conditions constantes qui ",
            "n'affectent pas significativement la probabilité de défaillance. Dans ce cas, la distribution de Weibull se réduit à une distribution exponentielle."
          ),
          strong("Pour plus d'informations au sujet de la distribution de Weibull:"),
          tags$ul(
            tags$li("Wikipedia : ", tags$a("Distribution de Weibull", href = "https://fr.wikipedia.org/wiki/Distribution_de_Weibull")),
            tags$li("Rinne, H. (2009). The Weibull Distribution: A Handbook. CRC Press."),
            tags$li("Forbes, C., Evans, M., Hastings, N., & Peacock, B. (2011). Statistical Distributions. John Wiley & Sons."),
            tags$li("Balakrishnan, N., & Basu, A. P. (2000). Statistical Methods in Reliability and Life Testing. World Scientific."),
            tags$li("Montgomery, J. C., & Runger, L. A. (1999). A Review of Weibull Statistics. Quality and Reliability Engineering International, 15(4), 275-291.")
          )
        ),
        
        
        tabPanel(
          strong(h4("Simulation")),
          fluidRow(
            column(
              width = 6,
              plotOutput("distPlot")
            ),
            column(
              width = 6,
              plotOutput("cdfPlot")
            )
          ),
          fluidRow(
            column(
              width = 6,
              conditionalPanel(
                condition = "input.input_mode == true",
                plotOutput("histPlot")
              ),
              conditionalPanel(
                condition = "input.input_mode == false",
                plotOutput("histPlot_sample")
              )
            ),
            column(
              width = 6,
              conditionalPanel(
                condition = "input.input_mode == true",
                plotOutput("cdfPlot_generated")
              ),
              conditionalPanel(
                condition = "input.input_mode == false",
                plotOutput("cdfPlot_sample")
              )
            )
          ),
          DT::DTOutput("statsOutput")
        ),
        tabPanel(
          strong(h4("Calcul")),
          fluidRow(
            column(
              width = 6,
              uiOutput("pdfEquation"),
              plotOutput("probabilityPlot")  # Add this line for the probability plot
            ),
            column(
              width = 6,
              uiOutput("cdfEquation"),
              plotOutput("quantilePlot") 
            )
          ),
          p(),
          p(),
          fluidRow(
            column(
              width = 4,
              numericInput("x_value", 
                           "Valeur de X:", 
                           value = 0, 
                           min = 0, 
                           step = 0.1,
                           width = "100px")
            ),
            column(
              width = 4, offset = 3,
              numericInput("prob_value", 
                           "Rang Centile:", 
                           value = 0.5, 
                           min = 0, 
                           max = 1, 
                           step = 0.01,
                           width = "100px")
            )
          ),
          fluidRow(
            column(
              width = 4,
              textOutput("computed_prob")
            ),
            column(
              width = 4, offset = 3,
              textOutput("computed_quantile")
            )
          ),
          p(),
          p(),
          fluidRow(
            column(
              width = 12,
              DTOutput("weibullChars")
            )
          )
          
        ),
        
        
        tabPanel(strong(h4("Sous R")),
                 p(),
                 
                 "Des fonctions R existent pour examiner différents aspects de la distribution de Weibull. En particulier, considérons les librairies ", strong("EnvStats"), " et ", strong("ForestFit"), ". ",
                 "La première offre les fonctions ", strong("eweibull"), " et ", strong("eqWeibull"), " pour estimer les paramètres de la distribution et ses quantiles, respectivement. ", p(),
                 "Pour illustrer, générons un ensemble de 500 observations tirées d'une distribution de Weibull dont les paramètres sont inconnus, et choisis au hasard: ", p(),

                 wellPanel(
                   style = "background-color: #f5f5f5; border: 1px solid #ddd; padding: 10px;",
                   tags$pre(
                     "# Données aléatoires, suivant un modèle
shape <- runif(1, 0.2, 10)
scale <- runif(1, 2, 500)
x <- rweibull(500, shape, scale)"
                   )
                 ),
                 
                 "À l'aide des fonctions disponibles, estimons les paramètres de la distribution d'où sont tirées les observations, en supposant qu'il s'agit bien d'une distribution de Weibull.", p(),

                 h4("EnvStats : eweibull et eqweibull"),
                 
                 h5(strong("Fonction eweibull")),
                 p("La fonction ", strong("EnvStats::eweibull"), " permet d'estimer les paramètres de la distribution de Weibull. "), 
                 p("Voici un exemple d'utilisation, utilisant les données préalablement générées.  On considère une distribution de Weibull à 2 paramètres, et on utilise une estimation par vraisemblance maximale:"),
 
                 wellPanel(
                   style = "background-color: #f5f5f5; border: 1px solid #ddd; padding: 10px;",
                   tags$pre(
                     "EnvStats::eweibull(x, method = 'mle')"
                   )
                 ),
                 
                 actionBttn("calc_eweibull", 
                            "Calculer les paramètres avec eweibull", 
                            color = "success", 
                            size = "sm",
                            style = "bordered"),
                 uiOutput("eweibull_result_container"),
                 
                 h5(strong("Fonction eqweibull")),
                 p("La fonction ", strong("EnvStats::eqweibull"), " permet non seulement d'estimer les paramètres de la distribution de Weibull, mais aussi d'en calculer les quantiles. On peut alimenter cette fonction ",
                 "àvec un objet produit par la fonction eweibull, ou un vecteur de données originales. Voici un exemple d'utilisation de cette fonction, en exploitant les résultats obtenus à l'aide de la fonction eweibull [parms]:"),
                 
                 wellPanel(
                   style = "background-color: #f5f5f5; border: 1px solid #ddd; padding: 10px;",
                   tags$pre(
                     "EnvStats::eqweibull(x = parms, 
                    p = c(0.10, 0.25, 0.5, 0.75, 0.90), 
                    method = 'mle', 
                    digits = 4)"
                   )
                 ),
                 
                 actionBttn("calc_eqweibull", 
                              "Calculer les quantiles avec eqweibull",
                              color = "success", 
                              size = "sm",
                              style = "bordered"),
                 uiOutput("eqweibull_result_container"),
                 
                 p(),
        p(),
        p(),
        h4("ForestFit : fitWeibull"),

                 p("La fonction ", strong("ForestFit::fitWeibull"), " permet d'examiner les distributions de Weibull impliquant deux [forme et échelle] ou trois  paramètres [forme, échelle et location/seuil]. Cette fonction nécessite que l'on propose des valeurs initiales pour chacun des paramètres [argument starts].", p(),
                 
                   h5(strong("Distribution de Weibull, 2 paramètres")),
                 "Pour les données générées précédemment, considérant une distribution de Weibull à deux paramètres, on propose forme = 1 et échelle = médiane/2. ",
                 "On obtient les résultats suivants: ", p(),
                 
                 wellPanel(
                   style = "background-color: #f5f5f5; border: 1px solid #ddd; padding: 10px;",
                   tags$pre(
                     "ForestFit::fitWeibull(x, 
                      location = FALSE, 
                      method = 'ml', 
                      starts = c(1, median(x)/2))"
                   )
                 ),

                 actionBttn("calc_fitWeibull2", 
                            "Calculer les paramètres avec fitWeibull [2 paramètres]",
                            color = "success", 
                            size = "sm",
                            style = "bordered"),
                 uiOutput("fitWeibull2_result_container"),
  
                 h5(strong("Distribution de Weibull, 3 paramètres")),
                 "Considérant une distribution de Weibull à trois paramètres, on propose les paramètres de forme et d'échelle obtenus précédemment, comme valeurs initiales, et on ajoute la valeur minimum de $x$ comme valeur initiale de la location [ou du seuil].",
                 "On obtient les résultats suivants: ", p(),
                   
                 wellPanel(
                   style = "background-color: #f5f5f5; border: 1px solid #ddd; padding: 10px;",
                   tags$pre(
          "parms2 <- fitWeibull(x, 
                     location = FALSE, 
                     method = 'ml', 
                     starts = c(1, median(x)/2))
                     
parms3 <- fitWeibull(x, 
                     location = TRUE,
                     method = 'mm1', 
                     starts = c(parms2$estimate, min(x)))"
                   )
                 ),
          
                   actionBttn("calc_fitWeibull3", 
                              "Calculer les paramètres avec fitWeibull [3 paramètres]",
                              color = "success", 
                              size = "sm",
                              style = "bordered"),
                   uiOutput("fitWeibull3_result_container")
                 ),
        
        h4("MASS : fitdistr"),
        
        p("La fonction ", strong("MASS::fitdistr"), " permet d'examiner les distributions de différents types, spécifiées par l'argument densefun. Cet argument peut prendre les valeurs suivantes: ",
        " beta, cauchy, chi-squared, exponential, gamma, geometric, lognormal, logistic, negative binomial, normal, Poisson, t et Weibull. ", p(),
        "Pour les données générées précédemment, on obtient: "),
          
        wellPanel(
          style = "background-color: #f5f5f5; border: 1px solid #ddd; padding: 10px;",
          tags$pre(
            "MASS::fitdistr(x, 
               densefun = 'weibull')
                     ")
        ),
        
        actionBttn("calc_fitdistr", 
                   "Calculer les paramètres avec fitdistr",
                   color = "success", 
                   size = "sm",
                   style = "bordered"),
        uiOutput("fitdistr_result_container")
        ),

        navbarMenu(h4("Exercices"),
                   tabPanel("Exercise 1",
                            h5("Exercise 1"),
                            p("La durée de vie des moteurs électriques suit une distribution de Weibull avec un paramètre de forme $\\alpha=1.5$ et d'échelle $\\beta=2000$ heures. Quelle est la probabilité qu'un moteur cesse de fonctionner avant la 1500ième heure d'utilisation?"),
                            actionBttn("btn_sol1", "Afficher la solution",
                                       color = "success", 
                                       size = "sm",
                                       style = "bordered"),
                            uiOutput("solution1")
                   ),
                   tabPanel("Exercise 2",
                            h5("Exercise 2"),
                            p("Dans les conditions décrite en [1], quelle est la probabilité qu'un moteur fonctionne encore après 4000 heures d'utilisation?"),
                            actionBttn("btn_sol2", "Afficher la solution",
                                       color = "success", 
                                       size = "sm",
                                       style = "bordered"),
                            uiOutput("solution2")
                   ),
                   tabPanel("Exercise 3",
                            h5("Exercise 3"),
                            p("Dans les conditions décrite en [1], quelle est la durée de vie que l'on peut s'attendre d'obtenir pour 90% des moteurs?"),
                            actionBttn("btn_sol3", "Afficher la solution",
                                       color = "success", 
                                       size = "sm",
                                       style = "bordered"),
                            uiOutput("solution3")
                   ),
                   tabPanel("Exercise 4",
                            h5("Exercise 4: Analyse Comparative"),
                            p(HTML("Deux types d'ampoules, A et B, sont envisagés pour une application particulière. L'ampoule A suit une distribution de Weibull avec des paramètres <code>(forme = 2, échelle = 1000 heures)</code>, tandis que l'ampoule B suit une distribution de Weibull avec des paramètres <code>(forme = 3, échelle = 1000 heures)</code>. Quel type d'ampoule est censé avoir un temps moyen de défaillance plus long ? Justifiez votre réponse.")),
                            actionBttn("btn_sol4", "Afficher la solution",
                                       color = "success", 
                                       size = "sm",
                                       style = "bordered"),
                            uiOutput("solution4")
                   ),

                   tabPanel("Exercise 5",
                            h5("Exercise 5: Estimation de paramètres"),
                            HTML(paste0("La commande <code>", strong("x <- rweibull(n = 5000, shape = runif(1, .5, 10), scale = runif(1, 100, 5000))"), "</code> génère 5000 observations suivant une distribution de Weibull dont les paramètres sont inconnus. Sous R, exécutez cette commande et estimez les paramètres de la distribution en utilisant l'une ou l'autre des fonctions proposées sous l'onglet ", strong("Sous R"), ".")),
                            p(),
                            actionBttn("btn_sol5", "Afficher la solution",
                                       color = "success", 
                                       size = "sm",
                                       style = "bordered"),
                            uiOutput("solution5")
                            )
                   )
        )
      )
    )
  )

server <- function(input, output, session) {
  observeEvent(input$reset_button, {
    updateSliderInput(session, "shape", value = 2)
    updateSliderInput(session, "scale", value = 1)
    
    updateNumericInput(session, "shape_num", value = 2)
    updateNumericInput(session, "scale_num", value = 1)
    
    updateSwitchInput(session, "input_type", value = TRUE)
    updateSwitchInput(session, "input_mode", value = TRUE)
    updateNumericInput(session, "num_observations", value = 100)
    
    output$histPlot_sample <- renderPlot(NULL)
    output$cdfPlot_sample <- renderPlot(NULL)
    output$cdfPlot_generated <- renderPlot(NULL)
    output$histPlot <- renderPlot(NULL)
    output$statsOutput <- DT::renderDT(NULL)
    output$fittedParams <- renderDT(NULL)
  })
  
  percentile_99 <- reactive({
    qweibull(0.99, shape = input$shape, scale = input$scale)
  })
  
  observe({
    max_x <- percentile_99()
    updateSliderInput(session, "x_value", max = max_x)
  })
  
  
  parms <- reactive({
    if (input$input_mode) {
      parms_estimated()
    } else {
      estimated_parms()
    }
    
  })
  
  estimated_parms <- reactiveVal(list(shape = 2, scale = 1)) 
  
  generate_observations <- eventReactive(input$generate_button, {
    params <- parms()
    observations <- rweibull(input$num_observations, shape = params$shape, scale = params$scale)
    return(observations)
  })
  
  data <- reactive({
    if (input$input_mode) {
      generate_observations()
    } else {
      req(input$file1)
      data <- read.csv(input$file1$datapath)
      if (!"x" %in% names(data) || !is.numeric(data$x)) {
        showNotification("The file must contain a numeric variable named 'x'.", type = "error")
        return(NULL)
      }
      return(data$x)
    }
  })
  
  observeEvent(input$file1, {
    req(input$file1)
    data <- read.csv(input$file1$datapath)
    if (!"x" %in% names(data) || !is.numeric(data$x)) {
      showNotification("The file must contain a numeric variable named 'x'.", type = "error")
      return(NULL)
    }
    fit <- fitdist(data$x, "weibull")
    estimated_parms(list(shape = fit$estimate[1], scale = fit$estimate[2]))
    
    updateSliderInput(session, "shape", value = fit$estimate[["shape"]])
    updateSliderInput(session, "scale", value = fit$estimate[["scale"]])
    
    updateNumericInput(session, "shape_num", value = fit$estimate[["shape"]])
    updateNumericInput(session, "scale_num", value = fit$estimate[["scale"]])
    
    output$fittedParams <- renderDT({
      formatted_params <- round(fit$estimate, 4)
      formatted_sds <- round(fit$sd, 4)
      datatable(
        data.frame(
          Parameter = c("Shape (\\alpha)", "Scale (\\beta)"),
          Estimate = formatted_params,
          Std.Err. = formatted_sds,
          stringsAsFactors = FALSE
        ),
        options = list(paging = FALSE, searching = FALSE, ordering = FALSE),
        rownames = FALSE
      )
    })
  })
  
  # Define a reactive expression for the estimated parameters
  parms_estimated <- reactive({
    if (input$input_type) {
      list(
        shape = input$shape,
        scale = input$scale
      )
    } else {
      list(
        shape = input$shape_num,
        scale = input$scale_num
      )
    }
  })
  
  output$pdfEquation <- renderText({
    params <- parms()
    paste("PDF: f(x) = (", params$shape, "/", params$scale, ") * (x/", params$scale, ")^(", params$shape - 1, ") * exp(-(x/", params$scale, ")^", params$shape, ")", sep = "")
  })
  
  output$cdfEquation <- renderText({
    params <- parms()
    paste("CDF: F(x) = 1 - exp(-(x/", params$scale, ")^", params$shape, ")", sep = "")
  })

  output$distPlot <- renderPlot({
    params <- parms()
    ggplot(data.frame(x = c(0, 100)), aes(x)) +
      stat_function(fun = dweibull, args = list(shape = params$shape, scale = params$scale)) +
      labs(
        title = "Weibull Distribution Density Function",
        x = "x",
        y = "Density"
      ) +
      theme_minimal() +
      scale_x_continuous(limits = c(0, qweibull(0.99, shape = params$shape, scale = params$scale)))
  })
  
  output$cdfPlot <- renderPlot({
    params <- parms()
    ggplot(data.frame(x = c(0, 10)), aes(x)) +
      stat_function(fun = pweibull, args = list(shape = params$shape, scale = params$scale)) +
      labs(
        title = "Weibull Distribution Cumulative Distribution Function",
        x = "x",
        y = "Cumulative Probability"
      ) +
      theme_minimal() +
      scale_x_continuous(limits = c(0, qweibull(0.99, shape = params$shape, scale = params$scale)))
  })
  
  output$histPlot <- renderPlot({
    req(data())
    bins <- length(hist(data(), breaks = "FD")$breaks)
    ggplot(data.frame(x = data()), aes(x)) +
      geom_histogram(aes(y = ..density..), bins = bins, fill = "lightblue", alpha = 0.5) +
      geom_density(color = "red", size = 0.7) +
      labs(
        title = "Histogram and Density Plot of Generated Data",
        x = "x",
        y = "Density"
      ) +
      theme_minimal()
  })
  
  output$cdfPlot_generated <- renderPlot({
    req(data())
    ecdf_data <- ecdf(data())
    ggplot(data.frame(x = c(c(0, qweibull(0.99, shape = input$shape, scale = input$scale)))), aes(x)) +
      stat_function(fun = function(x) ecdf_data(x)) +
      labs(
        title = "Cumulative Distribution of Generated Data",
        x = "x",
        y = "Cumulative Probability"
      ) +
      theme_minimal()
  })
  
  output$histPlot_sample <- renderPlot({
    req(data())
    bins <- length(hist(data(), breaks = "FD")$breaks)
    
    ggplot(data.frame(x = data()), aes(x)) +
      geom_histogram(aes(y = ..density..), bins = bins, fill = "lightblue", alpha = 0.5) +
      geom_density(color = "red", size = 0.8) +
      labs(
        title = "Distribution de l'échantillon",
        x = "x",
        y = "Densité"
      ) +
      theme_minimal()
  })
  
  
  output$cdfPlot_sample <- renderPlot({
    req(data())
    ecdf_data <- ecdf(data())
    ggplot(data.frame(x = c(0, qweibull(0.999, shape = input$shape, scale = input$scale))), aes(x)) +
      stat_function(fun = function(x) ecdf_data(x)) +
      labs(
        title = "Distribution Cumulative de l'échantillon",
        x = "x",
        y = "Probabilité cumulative"
      ) +
      theme_minimal()
  })
  
  output$statsOutput <- DT::renderDT({
    req(data())
    data_vals <- data()
    if (!is.null(data_vals)) {
      params <- parms()
      dist_params <- data.frame(
        Statistic = c("Forme (\u03B1)", "Échelle (\u03B2)"),
        Value = c(params$shape, params$scale)
      )
      sample_stats <- data.frame(
        Statistic = c("Moyenne", "Médiane", "Écart-Type", "Minimum", "Maximum"),
        Value = round(c(mean(data_vals), median(data_vals), sd(data_vals), min(data_vals), max(data_vals)), 4)
      )
      combined_stats <- rbind(dist_params, sample_stats)
      datatable(combined_stats, options = list(paging = FALSE, searching = FALSE, ordering = FALSE), rownames = FALSE)
    } else {
      NULL
    }
  })
  
  computed_pr <- reactive({
    params <- parms()
    x_value <- input$x_value
    prob <- pweibull(x_value, shape = params$shape, scale = params$scale)
    paste("P[X <", x_value, "] = ", round(prob, 4))
  })

  computed_quantile <- reactive({
    params <- parms()
    prob_value <- input$prob_value
    quantile <- qweibull(prob_value, shape = params$shape, scale = params$scale)
    paste("Q[", prob_value, "] = ", round(quantile, 4))
  })
  
  output$computed_prob <- renderText({
      computed_pr()
  })
  
  output$computed_quantile <- renderText({
    computed_quantile()
  })
  
  output$pdfEquation <- renderUI({
    params <- parms()
    withMathJax(
      helpText(
      "PDF Equation: $$f(x) = \\frac{\\alpha}{\\beta} \\left( \\frac{x}{\\beta} \\right)^{\\alpha-1} e^{-\\left( \\frac{x}{\\beta} \\right)^\\alpha}$$",
      "Where: $\\alpha$ =", round(params$shape, 6), ", $\\beta$ =", round(params$scale, 6)
      )  
    )
  })
  
  output$cdfEquation <- renderUI({
    params <- parms()
    helpText(
      "CDF Equation: $$F(x) = 1 - e^{-\\left( \\frac{x}{\\beta} \\right)^\\alpha}$$",
      "Where: $\\alpha$ =", round(params$shape, 6), ", $\\beta$ =", round(params$scale, 6)
    )
  })
  
  output$probabilityPlot <- renderPlot({
    params <- parms()
    data <- req(data()) 
    max_x <- max(data)
    
    x_values <- seq(0, max_x, by = 0.1)
    y_values <- dweibull(x_values, shape = params$shape, scale = params$scale)
    x_value <- input$x_value
    prob <- pweibull(x_value, shape = params$shape, scale = params$scale)
    
    df <- data.frame(x = x_values, y = y_values)
    shading_df <- data.frame(
      x = c(0, subset(df, x <= x_value)$x, x_value),
      y = c(0, subset(df, x <= x_value)$y, 0)
    )
    
    ggplot(df, aes(x = x, y = y)) +
      geom_line() +
      geom_polygon(data = shading_df, aes(x, y), fill = "lightblue", alpha = 0.5) +
      labs(x = "x", y = "Density", title = "Probability Density Function") +
      xlim(0, qweibull(0.99, shape = params$shape, scale = params$scale)) +
      geom_vline(xintercept = x_value, linetype = "dashed") +
      annotate("text", x = x_value, y = 0, vjust = 1, label = sprintf("x = %.2f", x_value), color = "red")
  })

  output$shapePlot <- renderPlot({
    shape_values <- c(0.5, 1, 1.5, 2, 3)
    x <- seq(0, 10, length.out = 1000)
    
    data <- data.frame(
      x = rep(x, times = length(shape_values)),
      shape = factor(rep(shape_values, each = length(x))),
      density = unlist(lapply(shape_values, function(k) dweibull(x, shape = k, scale = 1)))
    )
    
    ggplot(data, aes(x = x, y = density, color = shape)) +
      geom_line() +
      ylim(0, 1.5) +
      xlim(0, 4) +
      labs(title = "Variation de \u03B1",
           x = "x",
           y = "Densité",
           color = "Forme (\u03B1)") +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.5),
            legend.text = element_text(size = 12),   
            legend.title = element_text(size = 14),
            plot.title = element_text(size = 18))
  })
  

  
  output$quantilePlot <- renderPlot({
    params <- parms()
    prob_value <- input$prob_value
    
    if (is.na(prob_value) || prob_value < 0 || prob_value > 1) {
      return()
    }
    
    quantile_value <- qweibull(prob_value, shape = params$shape, scale = params$scale)
    
    prob_values <- seq(0, 1, by = 0.01)
    x_values <- qweibull(prob_values, shape = params$shape, scale = params$scale)
    
    valid_indices <- is.finite(x_values) & !is.na(x_values)
    x_values <- x_values[valid_indices]
    prob_values <- prob_values[valid_indices]
    
    if (length(x_values) == 0 || length(prob_values) == 0) {
      return()
    }
    
    df <- data.frame(x = x_values, prob = prob_values)
    
    ggplot(df, aes(x = x, y = prob)) +
      geom_line() +
      geom_hline(yintercept = prob_value, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = quantile_value, linetype = "dashed", color = "red") +
      geom_point(x = quantile_value, y = prob_value, color = "black", shape = 16) +
      labs(x = "Quantile", y = "Cumulative Probability", title = "Quantile Function") +
      ylim(0, 1) +
      xlim(range(x_values))
  })
  
  
  output$scalePlot <- renderPlot({
    scale_values <- c(0.5, 1, 1.5, 2, 3)
    x <- seq(0, 10, length.out = 1000)
    
    data <- data.frame(
      x = rep(x, times = length(scale_values)),
      scale = factor(rep(scale_values, each = length(x))),
      density = unlist(lapply(scale_values, function(lambda) dweibull(x, shape = 1.5, scale = lambda)))
    )
    
    ggplot(data, aes(x = x, y = density, color = scale)) +
      geom_line() +
      ylim(0, 1.5) +
      xlim(0, 10) +
      labs(title = "Variation de \u03B2",
           x = "x",
           y = "Densité",
           color = "Échelle (\u03B2)") +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.5),
            legend.text = element_text(size = 12),   
            legend.title = element_text(size = 14),
            plot.title = element_text(size = 18))
  })
  
  output$locationPlot <- renderPlot({
    location_values <- c(-1, 0, 1, 2)
    x <- seq(-2, 10, length.out = 1000)
    
    data <- data.frame(
      x = rep(x, times = length(location_values)),
      location = factor(rep(location_values, each = length(x))),
      density = unlist(lapply(location_values, function(mu) dweibull3(x, mu, beta = 1.2, lambda = 1.5)))
    )
    
    ggplot(data, aes(x = x, y = density, color = location)) +
      geom_line() +
      ylim(0, 0.6) +
      xlim(-2, 7.5) +
      labs(title = "Variation de \u03B3",
           x = "x",
           y = "Densité",
           color = "Location (\u03B3)") +
      theme_minimal() +
      theme(legend.position = c(0.85, 0.5),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            plot.title = element_text(size = 18))
  })

  output$weibullChars <- renderDT({
    library(latex2exp)
    shape <- input$shape
    scale <- input$scale
    characteristics <- weibull_characteristics(shape, scale)
    
    names(characteristics)[1] <- "Caractéristique"
    names(characteristics)[2] <- "Valeur"
    characteristics$Valeur <- round(characteristics$Valeur, 4)

    datatable(
      characteristics,
      options = list(
        pageLength = 12,
        autoWidth = TRUE,
        dom = 't' # Only show the table body
      ),
      rownames = FALSE,
      escape = FALSE # Allow HTML rendering of LaTeX equations
    )
  })
  
  # Create reactive values to keep track of the toggle state
  toggle_sol1 <- reactiveVal(FALSE)
  toggle_sol2 <- reactiveVal(FALSE)
  toggle_sol3 <- reactiveVal(FALSE)
  toggle_sol4 <- reactiveVal(FALSE)
  toggle_sol5 <- reactiveVal(FALSE)
  
  observeEvent(input$btn_sol1, {
    toggle_sol1(!toggle_sol1())
    if (toggle_sol1()) {
      updateActionButton(session, "btn_sol1", label = "Cacher la solution")
    } else {
      updateActionButton(session, "btn_sol1", label = "Afficher la solution")
    }
  })
  
  observeEvent(input$btn_sol2, {
    toggle_sol2(!toggle_sol2())
    if (toggle_sol2()) {
      updateActionButton(session, "btn_sol2", label = "Cacher la solution")
    } else {
      updateActionButton(session, "btn_sol2", label = "Afficher la solution")
    }
  })
  
  observeEvent(input$btn_sol3, {
    toggle_sol3(!toggle_sol3())
    if (toggle_sol3()) {
      updateActionButton(session, "btn_sol3", label = "Cacher la solution")
    } else {
      updateActionButton(session, "btn_sol3", label = "Afficher la solution")
    }
  })
  
  observeEvent(input$btn_sol4, {
    toggle_sol4(!toggle_sol4())
    if (toggle_sol4()) {
      updateActionButton(session, "btn_sol4", label = "Cacher la solution")
    } else {
      updateActionButton(session, "btn_sol4", label = "Afficher la solution")
    }
  })
  
  observeEvent(input$btn_sol5, {
    toggle_sol5(!toggle_sol5())
    if (toggle_sol5()) {
      updateActionButton(session, "btn_sol5", label = "Cacher la solution")
    } else {
      updateActionButton(session, "btn_sol5", label = "Afficher la solution")
    }
  })
  
  output$solution1 <- renderUI({
    if (toggle_sol1()) {
      withMathJax(
        helpText(
          strong("1. Méthode analytique"), ": ", p(),
          "$$F(x; \\alpha, \\beta) = 1 - e^{-\\left(\\frac{x}{\\beta}\\right)^\\alpha} = 1 - e^{-\\left(\\frac{1500}{2000}\\right)^{1.5}}=0.478 $$", p(),
          strong("2. Avec R"), ": ", p(),
          wellPanel(
            HTML("<pre><code>p <- pweibull(1500, shape = 1.5, scale = 2000)</code></pre>")
          ),
          strong("3. Avec le tableau de bord"), ": ", p(),
          wellPanel(
            "a. Dans le tableau de bord, spécifiez $\\alpha = 1.5$ et $\\beta = 2000$", br(),
            "b. Sous l'onglet ", strong("CALCUL"), "inscrivez la valeur de X = 1500", br(),
            "c. Prenez note du résultat..."
          ), p(),
          strong("4. Par simulation (R)"), " :", p(),
          wellPanel(
            HTML("<pre><code>x <- rweibull(n = 1000000, shape = 1.5, scale = 2000)
p <- mean(x <= 1500)
paste('p(x < 1500) = ', p)</code></pre>")
          )
        )
      )
    }
  })
  
  output$solution2 <- renderUI({
    if (toggle_sol2()) {
      withMathJax(
        helpText(
          strong("1. Méthode analytique"), ": ", p(),
          "$$F(x; \\alpha, \\beta) = 1 - e^{-\\left(\\frac{x}{\\beta}\\right)^\\alpha} = 1 - e^{-\\left(\\frac{4000}{2000}\\right)^{1.5}}=0.0591$$", p(),
          strong("2. Avec R"), ": ", p(),
          wellPanel(
            HTML("<pre><code>p <- pweibull(4000, shape = 1.5, scale = 2000, lower.tail = FALSE)</code></pre>")
          ),
          strong("3. Avec le tableau de bord"), ": ", p(),
          wellPanel(
            "a. Dans le tableau de bord, spécifiez $\\alpha = 1.5$ et $\\beta = 2000$", br(),
            "b. Sous l'onglet ", strong("CALCUL"), "inscrivez la valeur de X = 4000", br(),
            "c. Calculez 1 - le résultat affiché..."
          ), p(),
          strong("4. Par simulation (R)"), " :", p(),
          wellPanel(
            HTML("<pre><code>x <- rweibull(1000000, 1.5, 2000)
p <- mean(x >= 4000)
paste('p(x > 4000) = ', p)</code></pre>")
          )
        )
      )
    }
  })
  
  output$solution3 <- renderUI({
    if (toggle_sol3()) {
      withMathJax(
        helpText(
          strong("1. Méthode analytique"), ": ", p(),
          "$$Q(p; \\alpha, \\beta) = \\alpha \\left(-\\ln(1-p) \\right)^{\\frac{1}{\\beta}} = 1.5 \\left(-\\ln(1-0.90) \\right)^{\\frac{1}{2000}} = 3487.44$$"
          , p(),
          strong("2. Avec R"), ": ", p(),
          wellPanel(
            HTML("<pre><code>q <- qweibull(0.90, shape = 1.5, scale = 2000)</code></pre>")
          ),
          strong("3. Avec le tableau de bord"), ": ", p(),
          wellPanel(
            "a. Dans le tableau de bord, spécifiez $\\alpha = 1.5$ et $\\beta = 2000$", br(),
            "b. Sous l'onglet ", strong("CALCUL"), "inscrivez la valeur de Rang Centile = 0.90", br(),
            "c. Prenez note du résultat affiché..."
          ), p(),
          strong("4. Par simulation (R)"), " :", p(),
          wellPanel(
            HTML("<pre><code>x <- rweibull(1000000, 1.5, 2000)
q <- quantile(x, 0.90)
paste('Q(0.90) = ', q)</code></pre>")
          )
        )
      )
    }
  })
  
  output$solution4 <- renderUI({
    if (toggle_sol4()) {
      withMathJax(
        helpText(
          strong("1. Méthode analytique"), ": ", p(),
          "$$\\mu_A = \\beta \\cdot \\Gamma\\left(1 + \\frac{1}{\\alpha}\\right) = 1000 \\cdot \\Gamma\\left(1 + \\frac{1}{2}\\right)=886.227$$",
          "$$\\mu_B = \\beta \\cdot \\Gamma\\left(1 + \\frac{1}{\\alpha}\\right) = 1000 \\cdot \\Gamma\\left(1 + \\frac{1}{3}\\right)=892.980$$",
          print("Les ampoules B ont une durée de vie moyenne supérieure à celle des ampoules A.")
          , p(),
          strong("2. Avec R"), ": ", p(),
          wellPanel(
            HTML("<pre><code>weibull_mean <- function(shape, scale) {
scale * gamma(1 + 1 / shape)
}
AmpA <- weibull_mean(2, 1000)
AmpB <- weibull_mean(3, 1000)
ifelse(AmpB > AmpB, 'Les ampoules B ont une durée de vie plus longue!', 
    'Les ampoules B ont une durée de vie plus courte!')
    </code></pre>")
          ),
          strong("3. Avec le tableau de bord"), ": ", p(),
          wellPanel(
            "a. Dans le tableau de bord, spécifiez $\\alpha = 2$ et $\\beta = 1000$", br(),
            "b. Sous l'onglet ", strong("Simulation"), "prenez note de la moyenne de la distribution", br(),
            "c. Répétez a et b, en spécifiant $\\alpha = 3$ et prenez note de la moyenne de la distribution", br(),
            "d. La distribution dont la moyenne est supérieure est celle des ampoules dont la durée de vie moyenne est la plus élevée..."
          ), p(),
          strong("4. Par simulation (R)"), " :", p(),
          wellPanel(
            HTML("<pre><code>xA <- rweibull(1000000, 2, 1000)
xB <- rweibull(1000000, 3, 1000)
MxA <- mean(xA)
MxB <- mean(xB)
ifelse(MxB > MxA, 'Les ampoules B ont une durée de vie plus longue!', 
            'Les ampoules B ont une durée de vie plus courte!')
</code></pre>")
          )
        )
      )
    }
  })
  
  output$solution5 <- renderUI({
    if (toggle_sol5()) {
          helpText(
            strong("1. Génération des données"), ": ", p(),
            wellPanel(
              HTML("<pre><code>x <- rweibull(n = 5000, 
              shape = runif(1, .5, 10), 
              scale = runif(1, 100, 5000))</code></pre>")
            ),
            strong("2. Estimation des paramètres"), ": ", p(),
            wellPanel(
              "Plusieurs solutions sont possibles: ",
              HTML("<pre><code>
a. MASS::fitdistr(x, 'weibull'),
b. EnvStats::eqweibull(x, method = 'mle'),
c. ForestFit::fitWeibull(x, 
                         location = FALSE, 
                         method = 'ml', 
                         starts = c(1, 1))</code></pre>"),
            ), p()
          )
    }
    })
  
  # Generate random data once and store it in a reactive value
  random_data <- reactive({
    set.seed(1234)
    shape <- runif(1, 0.2, 10)
    scale <- runif(1, 2, 500)
    rweibull(500, shape, scale) 
  })
  
  rv1 <- reactiveValues(show_result = FALSE)
  rv2 <- reactiveValues(show_result = FALSE)
  
  observeEvent(input$calc_eweibull, {
    rv1$show_result <- !rv1$show_result  # Inverser l'état d'affichage
    
    if (rv1$show_result) {
      x <- random_data()
      parms <- eweibull(x, method = "mle")
      output$eweibull_result <- renderPrint({
        parms
      })
      updateActionButton(session, "calc_eweibull", label = "Cacher les résultats")
    } else {
      output$eweibull_result <- renderPrint({
        NULL
      })
      updateActionButton(session, "calc_eweibull", label = "Calculer les paramètres avec eweibull")
    }
  })
  
  # Output container for the result
  output$eweibull_result_container <- renderUI({
    if (rv1$show_result) {
      verbatimTextOutput("eweibull_result")
    } else {
      NULL  # Hide the output container when results are not shown
    }
  })
  
  observeEvent(input$calc_eqweibull, {
    rv2$show_result <- !rv2$show_result  # Inverser l'état d'affichage
    
    if (rv2$show_result) {
      x <- random_data()
      parms <- eweibull(x, method = "mle")
      qtls <- eqweibull(x = parms, 
                        p = c(0.10, 0.25, 0.5, 0.75, 0.90), 
                        method = "mle", 
                        digits = 4)
      output$eqweibull_result <- renderPrint({
        qtls
      })
      updateActionButton(session, "calc_eqweibull", label = "Cacher les résultats")
    } else {
      output$eqweibull_result <- renderPrint({
        ""  # Set to empty string to avoid showing NULL
      })
      updateActionButton(session, "calc_eqweibull", label = "Calculer les quantiles avec eqweibull")
    }
  })
  
  # Output container for the result
  output$eqweibull_result_container <- renderUI({
    if (rv2$show_result) {
      verbatimTextOutput("eqweibull_result")
    } else {
      NULL  # Hide the output container when results are not shown
    }
  })
  
  rv <- reactiveValues(show_result2 = FALSE, show_result3 = FALSE)
  
  # Observe button click event for 2 parameters
  observeEvent(input$calc_fitWeibull2, {
    rv$show_result2 <- !rv$show_result2  # Toggle display state
    
    if (rv$show_result2) {
      x <- random_data()
      parms2 <- fitWeibull(x, location = FALSE, method = "ml", starts = c(1, median(x)/2))
      output$fitWeibull2_result <- renderPrint({
        parms2
      })
      updateActionButton(session, "calc_fitWeibull2", label = "Cacher les résultats")
    } else {
      output$fitWeibull2_result <- renderPrint({
        NULL  # Ensure no display when hiding results
      })
      updateActionButton(session, "calc_fitWeibull2", label = "Calculer les paramètres avec fitWeibull (2 paramètres)")
    }
  })
  
  # Output container for the 2 parameters result
  output$fitWeibull2_result_container <- renderUI({
    if (rv$show_result2) {
      verbatimTextOutput("fitWeibull2_result")
    } else {
      NULL  # Hide the output container when results are not shown
    }
  })
  
  # Observe button click event for 3 parameters
  observeEvent(input$calc_fitWeibull3, {
    rv$show_result3 <- !rv$show_result3  # Toggle display state
    
    if (rv$show_result3) {
      x <- random_data()
      parms2 <- fitWeibull(x, location = FALSE, method = "ml", starts = c(1, median(x)/2))
      parms3 <- fitWeibull(x, location = TRUE, method = "mm1", starts = c(parms2$estimate, min(x)))
      output$fitWeibull3_result <- renderPrint({
        parms3
      })
      updateActionButton(session, "calc_fitWeibull3", label = "Cacher les résultats")
    } else {
      output$fitWeibull3_result <- renderPrint({
        NULL  # Ensure no display when hiding results
      })
      updateActionButton(session, "calc_fitWeibull3", label = "Calculer les paramètres avec fitWeibull (3 paramètres)")
    }
  })
  
  # Output container for the 3 parameters result
  output$fitWeibull3_result_container <- renderUI({
    if (rv$show_result3) {
      verbatimTextOutput("fitWeibull3_result")
    } else {
      NULL  # Hide the output container when results are not shown
    }
  })
  
  rvFitDistr <- reactiveValues(show_result = FALSE)
  
  observeEvent(input$calc_fitdistr, {
    rvFitDistr$show_result <- !rvFitDistr$show_result  # Inverser l'état d'affichage
    if (rvFitDistr$show_result) {
      x <- random_data()
      parms <- fitdistr(x, densfun = "weibull")
      output$FitDistr_result <- renderPrint({
      parms
      })
      updateActionButton(session, "calc_FitDistr", label = "Cacher les résultats")
    } else {
      output$FitDistr_result <- renderPrint({
        NULL
      })
      updateActionButton(session, "calc_FitDistr", label = "Calculer les paramètres avec fitdistr")
    }
    
    output$fitdistr_result_container <- renderUI({
      if (rvFitDistr$show_result) {
        verbatimTextOutput("FitDistr_result")
      } else {
        NULL
      }
    })
  })
  
  observeEvent(input$toggle_fitdistr, {
    toggle("fitdistr_result_container")
  })

  output$download_data <- downloadHandler(
    filename = function() {
      paste("weibull_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data.frame(x = data()), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
