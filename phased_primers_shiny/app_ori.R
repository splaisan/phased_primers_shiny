library(shiny)
### Define UI ----
ui <- fluidPage(
  titlePanel("Tool to design phased PCR1 amplicon primers"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("phasing", h6("define amount of primer phasing"),
                  min = 0, max = 21, value = 0),
      
      textInput("primer", h6("enter gene specific primer sequence (5'-> 3')"), 
                value = "CCTAHGGGRBGCAGCAG"),
      
      helpText(strong("standard primer sequences 16S V3/V4:"), 
               p("341-fw: CCTACGGGNGGCWGCAG"), p("805-re: GACTACHVGGGTATCTAATCC"),
               br()),
      textInput("adapter", h6("enter Adapter sequence sequence (5'-> 3')"), 
                value = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"),
      
      helpText(strong("adapter used at NGI:"), 
               p("fw: ACACTCTTTCCCTACACGACGCTCTTCCGATCT"), p("re: GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"),
               br()),
               
      textInput("prefix", h4("enter a prefix for primer names"),
                value = "primer"),
      br(),
               
      actionButton("update", h6("get primers")),
      downloadButton("downloadData", h6("export primers to csv file")),

      ),
    mainPanel(
      htmlOutput("infotext_upper"),
      imageOutput("prep_setup"),
      htmlOutput("infotext_lower"),
      plotOutput("plotting", inline = TRUE),
      tableOutput("primer_seq")

      
    )
  )
)
test <- "test"

### Define server logic ----
server <- function(input, output){
  adjust_special_nucleotides <- function(nucleotide_table){
    special_nucleotide_list <- compile_library_of_special_nucs()
    if(ncol(nucleotide_table) > 4){
      standard_nucleotides <- nucleotide_table[,c(1:4)]
      special_nucleotides <- as.data.frame(nucleotide_table[,c(5:ncol(nucleotide_table))])
      for (sp_nuc in special_nucleotides[2,]){
        sp_nuc <- as.character(sp_nuc)
        definition_sp_nuc <- special_nucleotide_list[sp_nuc][[1]]
       # definition_sp_nuc <- get_nuc_special_nucleotides(sp_nuc, special_nucleotide_list)
        for(nuc in definition_sp_nuc){
          occurance_sp_nuc <- as.numeric(as.character(as.data.frame(special_nucleotides[,special_nucleotides[2,]==sp_nuc])[1,]))
          fraction_to_add <- occurance_sp_nuc/length(definition_sp_nuc)
          standard_nucleotides[,standard_nucleotides[2,]==nuc][1] <- as.numeric(standard_nucleotides[,standard_nucleotides[2,]==nuc][1]) + fraction_to_add
        }
      }
      return(standard_nucleotides)
    }
    else{
      return(nucleotide_table)
    }
  }
  phase_primer <- function(primer, amount_of_phasing, special_nucleotide_library ){
    if(amount_of_phasing > 0){
      # I separate the primers by nucleotide
      split_primer <- strsplit(primer,"")
      # I copy the lists of nucleotides X times, where X is the number of additional phases
      lists_of_nucleotides <- split_primer
      for(i in 1:amount_of_phasing){
        lists_of_nucleotides <- c(lists_of_nucleotides, split_primer)
      }
      # Now I calculate the phasing nucleotides that will be added to the phases
      # I enter a loop that loops as many times as I have phases
      # for each loop occurance I count down
      phase_count <- amount_of_phasing
      nucleotide_occurance <- data.frame("A"=0,"T"=0,"C"=0,"G"=0)
      added_nucleotides <- c()
      for(i in 1:amount_of_phasing){
        # I take the subset of nucleotides from position 1 until the position of the count down
        subset_original_list_of_nucleotides <- lists_of_nucleotides[[1]][1:phase_count]
        # I calculate the occurcance of each nucleotide in the subset and put them into a df
        nucleotide_occurance_subset <- data.frame(table(c(as.character(added_nucleotides),subset_original_list_of_nucleotides)))
        nuc_occurance_df <- data.frame(t(nucleotide_occurance_subset[-1]))
        colnames(nuc_occurance_df) <- as.character(nucleotide_occurance_subset[,1])
        nucleotide_count_subset_with_special <- merge.data.frame(nucleotide_occurance, nuc_occurance_df, all.y=T, sort = F)
        # numbers for special nucleotides are added based on their proportion to the standard nucleotides
        # the nucleotide that occurs the least in that position is chosen
        nucleotide_count_subset_with_special[is.na(nucleotide_count_subset_with_special)] <- 0
        nucleotide_count_subset_with_special[2,] <- colnames(nucleotide_count_subset_with_special)
        nucleotide_count_subset <- adjust_special_nucleotides(nucleotide_count_subset_with_special)
        lowest_occurance <- as.data.frame(min(as.numeric(nucleotide_count_subset[1,])))
        nuc_lowest_occurance <- nucleotide_count_subset[2,nucleotide_count_subset[1,]==lowest_occurance[1,1]]
        while(length(nuc_lowest_occurance) != 1){
          # then up to X positions positions in row 1 (starting from 1,1) are added until the tie is broken
          # number of nucleotides to be excluded is one less as the length of the list of nucleotides with the lowest occurances
          low_occ_nuc_in_primer <- match(nuc_lowest_occurance, lists_of_nucleotides[[i]])
          low_occ_nuc_in_primer[is.na(low_occ_nuc_in_primer)]<- 0
          if(min(low_occ_nuc_in_primer) < length(subset_original_list_of_nucleotides)){   
            position_to_remove <- which(low_occ_nuc_in_primer==min(low_occ_nuc_in_primer))
            if(length(position_to_remove) != length(low_occ_nuc_in_primer)){
              nuc_lowest_occurance <- nuc_lowest_occurance[- position_to_remove]
            }
            else{
              nuc_lowest_occurance <- sample(nuc_lowest_occurance, 1)
            }
          }
          else{
            nuc_lowest_occurance <- sample(nuc_lowest_occurance, 1)
          }
        }
        # the chosen nucleotide is remembered
        added_nucleotides <- c(nuc_lowest_occurance, added_nucleotides)
        # the chosen nucleotide is added in front of the list that records phasing nucleotides
        rows_with_nuc_to_add <- c((amount_of_phasing-phase_count+2):(amount_of_phasing+1))
        for(row_to_add_nuc in rows_with_nuc_to_add){
          lists_of_nucleotides[[row_to_add_nuc]] <- c(as.character(nuc_lowest_occurance),lists_of_nucleotides[[row_to_add_nuc]])
        }
        # The count down number is reduced by one
        phase_count <- phase_count - 1  
      } 
    }
    else{
    lists_of_nucleotides <- strsplit(primer,"")
    }
    return(lists_of_nucleotides)
  }
  design_primers <- function(phased_primer, adapter, amount_of_phasing){
    full_phased_primer <- unlist(lapply(phased_primer(), function(x){ paste(c(adapter,x), collapse = '')}))
    primer_names <- paste(input$prefix, c(1:length(full_phased_primer)), sep="_")
    df_primers <- data.frame(primer_name = primer_names, sequence = full_phased_primer)
    return(df_primers)
  }
  compile_library_of_special_nucs <- function(){
    special_nucleotides_list <- list()
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list, "R", c("A","G"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"Y" , c("C", "T"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"S" , c ("G", "C"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"W" , c("A","T"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"K" , c("G", "T"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"M" , c("A", "C"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"B" , c("C", "G", "T"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"D" , c("A", "G", "T"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"H" , c("A", "C", "T"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"V" , c("A", "C", "G"))
    special_nucleotides_list <-add_lis_to_listolists(special_nucleotides_list,"N" , c("A", "C", "G", "T"))
    return(special_nucleotides_list)
  }
  add_lis_to_listolists <- function(list_of_lists, new_list_name, new_list_content){
    list_of_lists[[new_list_name]] <- new_list_content
    return(list_of_lists)
  }
  plot_nucleotide_abundance_phased <- function(list_phased_primer, length_to_plot, title){
    short_phased_primers <- unlist(lapply(list_phased_primer, function(X){return(X[1:length_to_plot])}))
    short_phased_split_primers <- matrix(short_phased_primers,nrow = length(list_phased_primer), ncol = length_to_plot, byrow=T)
    rel_abundance_by_position <- list()
    for (i in c(1:ncol(short_phased_split_primers))){
      table_nuc_list <- table(short_phased_split_primers[,i])
      rel_abundance_list <- calculate_relative_nuc_abundance(table_nuc_list)[1,]
      perc_abundance_list <- as.data.frame(100/sum(as.numeric(rel_abundance_list[1,]))*as.numeric(rel_abundance_list[1,]), row.names=colnames(rel_abundance_list))
      rel_abundance_by_position <- add_lis_to_listolists(rel_abundance_by_position,i, t(perc_abundance_list))
    }
    matrix_rel_abundance_by_position <- do.call(rbind, lapply(rel_abundance_by_position, function(x) x[match(colnames(rel_abundance_by_position[[1]]), colnames(x))]))
    colnames(matrix_rel_abundance_by_position) <- colnames(rel_abundance_by_position[[1]])
    class(matrix_rel_abundance_by_position) <- "numeric"
    rownames(matrix_rel_abundance_by_position) <- as.character(c(1:12))
    bp_primer <- barplot(t(matrix_rel_abundance_by_position), beside = TRUE, legend.text = T,col=c("forestgreen","red", "cornflowerblue", "black"), 
                         space = c(0,2),args.legend = list(x = "topright", bty = "n", inset=c(-0.06, 0)), ylim =c(0,100), 
                         xlab= "position", ylab="relative abundance in %", cex.names = 0.75, main=title)
    return(bp_primer)
  }
  calculate_relative_nuc_abundance <- function(table_nuc_list, special_nucleotides_list){
    nucleotide_occurance <- data.frame("A"=0,"T"=0,"C"=0,"G"=0)
    
    nucleotide_occurance_subset <- data.frame(table_nuc_list)
    nuc_occurance_df <- data.frame(t(nucleotide_occurance_subset[-1]))
    colnames(nuc_occurance_df) <- as.character(nucleotide_occurance_subset[,1])
    nucleotide_count_subset_with_special <- merge.data.frame(nucleotide_occurance, nuc_occurance_df, all.y=T, sort = F)
    # numbers for special nucleotides are added based on their proportion to the standard nucleotides
    # the nucleotide that occurs the least in that position is chosen
    nucleotide_count_subset_with_special[is.na(nucleotide_count_subset_with_special)] <- 0
    nucleotide_count_subset_with_special[2,] <- colnames(nucleotide_count_subset_with_special)
    nucleotide_count_subset <- adjust_special_nucleotides(nucleotide_count_subset_with_special)
    
    return(nucleotide_count_subset)
  }
  
  ##Delay
  primer_delay <- eventReactive(input$update, {input$primer})
  adapter_delay <- eventReactive(input$update, {input$adapter})
  prefix_delay <- eventReactive(input$update, {input$prefix})
  observeEvent(input$update, {infotext_upper(input$update)})
  observeEvent(input$update, {infotext_lower("")})
  
  observeEvent(input$update, {image_prep_setup(list(src = "www/alt.png", width = 20, height =20))})
 # observeEvent(input$update,{})
  ##End delay
  
  #Define run
  infotext_upper <- reactiveVal(
    value = 
      "<h3><b>test<b></h3>
        <p>test2</p>")
  infotext_lower <- reactiveVal(
    value = 
      "<h3><b>test3<b></h3>
    <p>test4</p>")
  image_prep_setup <-   reactiveVal( value = list(src = "www/amplicon_seq.png", width=40, height = 50))
  
  designed_phase <- reactive({phase_primer(primer_delay(), input$phasing, special_nucleotides_list)})
  assembled_primers <- reactive({design_primers(designed_phase, adapter_delay(), input$phasing)})
  plot_primer <- reactive({plot_nucleotide_abundance_phased(designed_phase(), 12, "nucleotide balance by position of sequencing")})
  
  ##Define output
  #output$infotext_upper <- renderText({infotext_upper()})
  #output$infotext_lower <- renderText({infotext_upper()})
  
  output$primer_seq <- renderTable({assembled_primers()})
  output$plotting <- renderPlot({plot_primer()},height = 400, width = 600)
  output$prep_setup <- renderImage({image_prep_setup()}, deleteFile = FALSE)
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      "phased_pcr-primers.csv"
    },
    content = function(file) {
      write.csv(assembled_primers(), file, row.names = FALSE)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
