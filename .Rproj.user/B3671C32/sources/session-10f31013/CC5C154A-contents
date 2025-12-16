library(shiny)
library(bslib)
library(shinyWidgets)
library(shiny.fluent)
library(shinyjs) 
library(shinycssloaders)
library(NGLVieweR)
library(DT)
library(dplyr)
###set AbRSA and tmp dir
AbRSA_path <- "~/Antigen_design/scripts/AbRSA/AbRSA"
tmp_dir <- "/data/sda/wt/tmp"
##
run_anno <- function(seq){
  temp_file <- paste0(tmp_dir,"/anno_test",".fa")
  seqinr::write.fasta(seq,names = "test",file.out = temp_file)
  commd <- paste0(AbRSA_path," -i ", temp_file, 
                  " -o ", tmp_dir, "/anno_test",
                  " > ", tmp_dir, "/anno_test",".out")
  system(commd)
  out <- tryCatch({
    data.table::fread(paste0(tmp_dir,"/anno_test",".out"),
                      skip = "#similarity",
                      data.table = F,sep = ":",header = F)
  },error = function(e){
    NA
  })
  return(out)
}

check_ab <- function(aaseq){
  ##输入的是H L A
  hla <- gsub(">[A-Z]","",aaseq) %>% 
    paste0("\n",.) %>% strsplit(.,"\n\n") %>% 
    `[[`(1) %>% `[`(2:4)
  hano <- run_anno(hla[1])
  lano <- run_anno(hla[2])
  if ((nrow(hano) >= 7) & (nrow(lano) >= 7)){
    return("yes")
  }else{
    return("no")
  }
}

check_abag_chains <- function(pdb_file, chain){
  pdbs <- bio3d::read.pdb(pdb_file)
  all_chains <- unique(pdbs$atom$chain)
  return(chain %in% all_chains)
}

link_docker <- tags$a(
  shiny::icon("docker"), "Docker",
  href = "https://github.com/rstudio/shiny",
  target = "_blank"
)
init_fasta <- ">H
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
>L
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
>A
CHPECFGPEADQCVACAHYKDPPFCVARCPSIWKFPDEEGACQPCPIN"

ui <- page_navbar(
  title = "Ab-Design",
  nav_panel(title = "Sequence",
            shinyFeedback::useShinyFeedback(),
            page_sidebar(
              sidebar = sidebar(
                title = "Settings",
                width = 300,
                selectInput("method","Design method",c("IgLM","AbGPT","ReprogBERT","PALM")),
                numericInput("design_num", "Number of Design Sequences", 30),
                selectInput("pre_structure","Antibody Structure Prediction",
                            c("TRUE","FALSE")),
                selectInput("docking","Antibody-Antigen Docking",
                            c("TRUE","FALSE")),
                PrimaryButton.shinyInput("run_seq","Run")
              ),
              card(
                height = "100px",
                card_header("Input"),
                fluidRow(
                  column(6, 
                         fluidRow(
                           textAreaInput("in_fasta",
                                         "Fasta Input (Three sequences, IDs are H L A)",
                                         width = "100%", height = "150px")
                         ),
                         fluidRow(
                           column(4,
                                  PrimaryButton.shinyInput(
                                    "load_example_fasta",
                                    text = "Load Example"
                                  )),
                           column(6,
                                  PrimaryButton.shinyInput(
                                    "reset_seq",
                                    text = "Reset Input"
                                  ))
                         )
                        ),
                  column(6, uiOutput("ag_pdb_upload"))
                )
              ),
              card(
                height = "100px",
                card_header("Output"),
                fluidRow(
                  column(4, textOutput("seq_pro_name")),
                  column(4, uiOutput("pre_struc_method")), ###展示哪种方法的
                  column(4, uiOutput("pre_struc_id")) ###展示哪个 ID 的
                  ),
                ##项目信息，展示哪些结构
                fluidRow(
                  column(6, DTOutput("seq_res",width = "500px"))
                ), ##表格，和展示PDB
                fluidRow(
                  column(8, NGLVieweROutput("seq_pre_struc_pdb")),
                  column(4, uiOutput("seq_down_ui")) ##下载结果
                )
              )
            )
          ),
  nav_panel(title = "Structure", 
            shinyFeedback::useShinyFeedback(),
            page_sidebar(
              sidebar = sidebar(
                title = "Settings",
                width = 300,
                selectInput("struc_method","Design method",
                            c("DiffAb","AbOpt","AbDockgen","RFantibody")),
                numericInput("struc_design_num", "Number of Design Sequences", 5),
                selectInput("pre_affinity","Antibody-Antigen Binding Affinity",
                            c("TRUE","FALSE")),
                selectInput("pre_energy","Antibody-Antigen Rosetta Energy",
                            c("TRUE","FALSE")),
                PrimaryButton.shinyInput("run_struc","Run")
              ),
              card(
                height = "100px",
                card_header("Input"),
                fluidRow(
                  column(6, 
                         fluidRow(
                           uiOutput("abag_pdb_upload")
                         ),
                         fluidRow(
                           column(4,
                                  PrimaryButton.shinyInput(
                                    "load_example_struc",
                                    text = "Load Example"
                                  )),
                           column(6,
                                  PrimaryButton.shinyInput(
                                    "reset_struc",
                                    text = "Reset Input"
                                  ))
                         )
                  ),
                  column(6, 
                         column(6,textInput("h_chain","The Heavy chain ID")),
                         column(6,textInput("ag_chain","The Antigen chain ID")))
                )
              ),
              card(
                height = "100px",
                card_header("Output"),
                fluidRow(
                  column(4, textOutput("struc_pro_name")),
                  column(4, uiOutput("struc_id_show")) ###展示哪个 ID 的
                ),
                ##项目信息，展示哪些结构
                fluidRow(
                  column(6, DTOutput("struc_res",width = "500px"))
                ), ##表格，和展示PDB
                fluidRow(
                  column(8, NGLVieweROutput("struc_design_pdb")),
                  column(4, uiOutput("struc_down_ui")) ##下载结果
                )
              )
            )
          ),
  nav_spacer(),
  nav_menu(
    title = "Use Docker",
    align = "right",
    nav_item(link_docker)
  ),
  theme = bs_theme(version = 5, bootswatch = "cerulean")
)

server <- function(input, output, session) {
  sc_path <- dirname(this.path::this.path())
  
  input_seq_data <- reactiveValues(
    fasta_path = NULL,
    ag_pab_path = NULL,
    res_dt = NULL,
    pro_name = NULL,
    ag_chain = NULL
  )
  output$ag_pdb_upload <- renderUI({
    fileInput("ag_pdb",
              "Antigen PDB File",
              placeholder = "No file selected.")
  })
  observeEvent(input$load_example_fasta,{
    output$ag_pdb_upload <- renderUI({
      fileInput("ag_pdb",
                "Antigen PDB File",
                placeholder = "Example File Ready.")
    })
    updateTextAreaInput(session, "in_fasta",value = init_fasta)
    input_seq_data$fasta_path <- "/home/ab/run/run_test/her2_abag.fasta"
    input_seq_data$ag_pab_path <- "/home/ab/run/run_test/her2_ag.pdb"
    input_seq_data$ag_chain <- "C"
  })
  observeEvent(input$reset_seq,{
    output$ag_pdb_upload <- renderUI({
      fileInput("ag_pdb",
                "Antigen PDB File",
                placeholder = "No file selected.")
    })
    updateTextAreaInput(session, "in_fasta",value = "")
    input_seq_data$fasta_path <- NULL
    input_seq_data$ag_pab_path <- NULL
  })
  ###run
  observeEvent(input$run_seq, {
    ###清空output
    system("docker exec --user ab abdesigner_last /bin/bash -c 'rm -rf /home/ab/run/run_test/output/*'")
    system(paste0("rm -rf ", sc_path, "/test/output/*"))
    ###
    tmp_file_pre <- as.character(as.numeric(Sys.time())) %>% gsub("[.]+","",.)
    input_seq_data$pro_name <- tmp_file_pre
    # 验证输入是否为H,L链
    has_seq_fasta <- check_ab(input$in_fasta) == "yes"
    shinyFeedback::feedbackDanger("in_fasta", !has_seq_fasta, "Input is not valid H or L")
    req(has_seq_fasta, cancelOutput = TRUE)
    ##保存成fasta文件
    hla <- gsub(">[A-Z]","",input$in_fasta) %>% 
      paste0("\n",.) %>% strsplit(.,"\n\n") %>% 
      `[[`(1) %>% `[`(2:4)
    if (is.null(input_seq_data$fasta_path)){
      fasta_name <- paste0("tmp_",tmp_file_pre,".fasta")
      input_seq_data$fasta_path <- paste0("/home/ab/run/run_test/output/",fasta_name)
      seqinr::write.fasta(as.list(hla),c("H","L","A"),
                          paste0(sc_path,"/test/output/",fasta_name))
    }
    if (is.null(input_seq_data$ag_pab_path)){
      agpdb_name <- paste0("tmp_",tmp_file_pre,".pdb")
      fs::file_copy(input$ag_pdb$datapath, paste0(sc_path,"/test/output/",agpdb_name))
      input_seq_data$ag_pab_path <- paste0("/home/ab/run/run_test/output/",agpdb_name)
      tmp_pdb <- bio3d::read.pdb(input$ag_pdb$datapath)
      input_seq_data$ag_chain <- unique(tmp_pdb$atom$chain)
    }
    ##生成配置文件
    config_seq <- tomledit::toml(
      input = list(
        fasta_file = input_seq_data$fasta_path,
        ag_pdb_file = input_seq_data$ag_pab_path
      ),
      settings = list(
        design_method = input$method,
        design_num = input$design_num,
        pre_structure = as.logical(input$pre_structure),
        docking = as.logical(input$docking),
        out_path = "/home/ab/run/run_test/output/"
      )
    )
    tomledit::write_toml(config_seq,
                         paste0(sc_path,"/test/output/tmp_seq_",tmp_file_pre,".toml"))
    ###run
    showPageSpinner()
    system(paste0("docker exec --user ab abdesigner_last /bin/bash -c 'source /home/ab/miniconda3/bin/activate; python /home/ab/run/run_seq.py -c ",
                  "/home/ab/run/run_test/output/tmp_seq_",tmp_file_pre,".toml'"))
    hidePageSpinner()
    ##打包结果
    zip::zipr(zipfile = paste0(sc_path,"/test/res/",tmp_file_pre,"_res.zip"), 
              files = paste0(sc_path,"/test/output/"))
    ###展示结果
    seq_dt <- read.csv(paste0(sc_path, "/test/output/seq_model_res.csv"))
    input_seq_data$res_dt <- seq_dt
    output$seq_pro_name <- renderText({
      paste0("Project Name:\n",tmp_file_pre)
    })
    output$seq_res <- renderDT(
      if (as.logical(input$pre_structure)){
        if (as.logical(input$docking)){
          seq_dt %>% dplyr::select(-X) %>% dplyr::select(-combined_id) %>% 
            dplyr::select(ID,Method,lDDT,pTM,HCDR3_RMSD,Docking_score,Sequence)
        }else{
          seq_dt %>% dplyr::select(-X) %>% dplyr::select(-combined_id) %>% 
            dplyr::select(ID,Method,lDDT,pTM,HCDR3_RMSD,Sequence)
        }
      }else{
        seq_dt %>% dplyr::select(-X) %>% dplyr::select(-combined_id) %>% 
          dplyr::select(ID,Method,Sequence)
      }, 
      options = list(pageLength = 5, dom = 'ft', searching = FALSE)
    )
    if (as.logical(input$pre_structure)){
      if (as.logical(input$docking)){
        output$pre_struc_method <- renderUI({
          selectInput("pre_struc_method_which","Which method to show?",
                      c("lgFold","tFold","Docking"),selected = "lgFold")
        })
        if (!dir.exists(paste0(sc_path, "/test/output/docking/"))){
          system(paste0("mkdir ",sc_path,"/test/output/docking;",
                        "unzip ",sc_path,"/test/output/docking_res",
                        ".zip -d ",sc_path,"/test/output/docking"))
        }
      }else{
        output$pre_struc_method <- renderUI({
          selectInput("pre_struc_method_which","Which method to show?",
                      c("lgFold","tFold"),selected = "lgFold",width="100%")
        })
      }
      
      output$pre_struc_id <- renderUI({
        selectInput("pre_struc_id_which","Which Seq ID to show?",
                    input_seq_data$res_dt$combined_id)
      })
      
      if (!dir.exists(paste0(sc_path,"/test/output/lgfold"))){
        system(paste0("mkdir ",sc_path,"/test/output/lgfold;",
                      "unzip ",sc_path,"/test/output/lgfold_pdb_*",
                      ".zip -d ",sc_path,"/test/output/lgfold"))
      }
      if (!dir.exists(paste0(sc_path,"/test/output/tfold"))){
        system(paste0("mkdir ",sc_path,"/test/output/tfold;",
                      "unzip ",sc_path,"/test/output/tfold_pdb_*",
                      ".zip -d ",sc_path,"/test/output/tfold"))
      }
    }
    output$seq_down_ui <- renderUI({
      downloadButton('seq_res_download', 'Download')
    })
  })
  output$seq_pre_struc_pdb <- renderNGLVieweR({
    req(input_seq_data$res_dt)
    req(as.logical(input$pre_structure))
    if (input$pre_struc_method_which == "lgFold"){
      pdb_path <- paste0(sc_path,"/test/output/lgfold/",
                         input$pre_struc_id_which,".pdb")
    }
    if (input$pre_struc_method_which == "tFold"){
      pdb_path <- paste0(sc_path,"/test/output/tfold/",
                         input$pre_struc_id_which,".pdb")
    }
    if (input$pre_struc_method_which == "Docking"){
      pdb_path <- paste0(sc_path,"/test/output/docking/",
                         input$pre_struc_id_which,"/",input$pre_struc_id_which,
                         "_rank_1_unrelaxed.pdb")
    }
    NGLVieweR(pdb_path) %>% 
      stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
      addRepresentation("cartoon",
                        param = list(name = "cartoon", colorScheme = "chainname")
      ) %>% 
      addRepresentation("surface", param = list(sele = paste0(":",input_seq_data$ag_chain)))
  })

  output$seq_res_download <- downloadHandler(
    filename = function() {
      paste0(input_seq_data$pro_name,"_res.zip")
    },
    content = function(file) {
      zip_file <- paste0(sc_path,"/test/res/",input_seq_data$pro_name,"_res.zip")
      file.copy(zip_file, file)
    }
  )
  #####structure
  input_struc_data <- reactiveValues(
    abag_pab_path = NULL,
    res_dt = NULL,
    pro_name = NULL
  )
  output$abag_pdb_upload <- renderUI({
    fileInput("abag_pdb",
              "Antibody-Antigen Complex PDB File",
              placeholder = "No file selected.")
  })
  observeEvent(input$load_example_struc,{
    output$abag_pdb_upload <- renderUI({
      fileInput("abag_pdb",
                "Antibody-Antigen Complex PDB File",
                placeholder = "Example File Ready.")
    })
    input_struc_data$abag_pab_path <- "/home/ab/run/run_test/output/1n8z_B_C.pdb"
    updateTextInput(inputId = "h_chain",value = "B")
    updateTextInput(inputId = "ag_chain",value = "C")
  })
  observeEvent(input$reset_struc,{
    output$abag_pdb_upload <- renderUI({
      fileInput("abag_pdb",
                "Antibody-Antigen Complex PDB File",
                placeholder = "No file selected.")
    })
    input_struc_data$abag_pab_path <- NULL
    updateTextInput(inputId = "h_chain",value = "")
    updateTextInput(inputId = "ag_chain",value = "")
  })
  ############
  observeEvent(input$run_struc,{
    ###清空output
    system("docker exec --user ab abdesigner_last /bin/bash -c 'rm -rf /home/ab/run/run_test/output/*'")
    system(paste0("rm -rf ",sc_path,"/test/output/*"))
    ##
    if (!is.null(input_struc_data$abag_pab_path)){
      if (input_struc_data$abag_pab_path == "/home/ab/run/run_test/output/1n8z_B_C.pdb"){
        system(paste0("cp ",sc_path,"/test/1n8z_B_C.pdb ",sc_path,"/test/output/1n8z_B_C.pdb")) 
      }
    }
    tmp_file_pre <- as.character(as.numeric(Sys.time())) %>% gsub("[.]+","",.)
    input_struc_data$pro_name <- tmp_file_pre
    if (is.null(input_struc_data$abag_pab_path)){
      abagpdb_name <- paste0("tmp_",tmp_file_pre,".pdb")
      fs::file_copy(input$abag_pdb$datapath, 
                    paste0(sc_path,"/test/output/",abagpdb_name))
      input_struc_data$abag_pab_path <- paste0("/home/ab/run/run_test/output/",abagpdb_name)
    }
    req(input_struc_data$abag_pab_path, cancelOutput = TRUE)
    ###
    own_path <- gsub("/home/ab/run/run_test/",
                     paste0(sc_path,"/test/"),input_struc_data$abag_pab_path)
    print(own_path)
    right_h <- check_abag_chains(own_path,input$h_chain)
    shinyFeedback::feedbackDanger("h_chain", !right_h, "H chain ID is not in PDB")
    right_ag <- check_abag_chains(own_path,input$ag_chain)
    shinyFeedback::feedbackDanger("ag_chain", !right_ag, "Antigen chain ID is not in PDB")
    req(right_ag, cancelOutput = TRUE)
    req(right_h, cancelOutput = TRUE)
    ##生成配置文件
    config_struc <- tomledit::toml(
      input = list(
        pdb_file = input_struc_data$abag_pab_path
      ),
      settings = list(
        design_method = input$struc_method,
        design_num = input$struc_design_num,
        Heavy_Chain = input$h_chain,
        Antigen_Chain = input$ag_chain,
        binding_affinity = as.logical(input$pre_affinity),
        binding_energy = as.logical(input$pre_energy),
        out_path = "/home/ab/run/run_test/output/"
      )
    )
    tomledit::write_toml(config_struc,
                         paste0(sc_path,"/test/output/tmp_struc_",tmp_file_pre,".toml"))
    ###
    ###run
    showPageSpinner()
    system(paste0("docker exec --user ab abdesigner_last /bin/bash -c 'source /home/ab/miniconda3/bin/activate; cd /home/ab/run/; python run_struc.py -c ",
                  "/home/ab/run/run_test/output/tmp_struc_",tmp_file_pre,".toml'"))
    hidePageSpinner()
    ##打包结果
    zip::zipr(zipfile = paste0(sc_path,"/test/res/",tmp_file_pre,"_res.zip"),
              files = paste0(sc_path,"/test/output/"))
    ###展示结果
    ##解压结果
    system(paste0("mkdir ",sc_path,"/test/output/res/;",
                  "unzip ",sc_path,"/test/output/*.zip",
                  " -d ",sc_path,"/test/output/res"))
    ##
    struc_dt <- read.csv(paste0(sc_path,"/test/output/struc_gen_seq.csv")) %>% select(-X)
    if (input$struc_method %in% c("DiffAb","AbOpt")){
      struc_dt$ID <- sprintf("%04d", struc_dt$ID)
      if (input$struc_method == "AbOpt"){
        struc_dt$ID <- paste0("gen_",struc_dt$ID)
      }
    } 
    if (input$struc_method == "RFantibody"){
      struc_dt$ID <- paste0(struc_dt$ID,"_dldesign_0")
      rf_score <- read.csv(paste0(sc_path,"/test/output/res/rfantibody_scores.csv")) %>% 
        select(-X)
      rf_score$ID <- paste0(rf_score$ID,"_dldesign_0")
      struc_dt <- left_join(struc_dt, rf_score)
    } 
    ##读取其他文件
    if (as.logical(input$pre_affinity) & (input$struc_method != "RFantibody")){
      aff <- read.csv(list.files(paste0(sc_path,"/test/output/res/"),
                                 pattern = "affinity",full.names = T)) %>% select(-X)
      if (input$struc_method == "DiffAb"){
        aff$id <- sprintf("%04d", aff$id)
      } 
      struc_dt <- left_join(struc_dt, aff %>% rename(ID=id))
    }
    if (as.logical(input$pre_energy)){
      energy <- read.csv(list.files(paste0(sc_path,"/test/output/res/"),
                                 pattern = "rosetta",full.names = T)) %>% select(-X)
      if (input$struc_method == "DiffAb"){
        energy$id <- sprintf("%04d", energy$id)
      } 
      struc_dt <- left_join(struc_dt, energy %>% rename(ID=id))
    }
    input_struc_data$res_dt <- struc_dt
    output$struc_pro_name <- renderText({
      paste0("Project Name:\n",tmp_file_pre)
    })
    output$struc_id_show <- renderUI({
      selectInput("struc_id_which","Which ID to show?",
                  input_struc_data$res_dt$ID)
    })
    output$struc_res <- renderDT(
      input_struc_data$res_dt, 
      options = list(pageLength = 5, dom = 'ft', searching = FALSE)
    )
    ##
    output$struc_down_ui <- renderUI({
      downloadButton('struc_res_download', 'Download')
    })
  })
  ###show pdb
  output$struc_design_pdb <- renderNGLVieweR({
    req(input_struc_data$res_dt)
    pdb_path <- paste0(sc_path,"/test/output/res/",
                       input$struc_id_which,".pdb")
    show_agchain <- case_when(
      input$struc_method == "AbDockgen" ~ "A",
      input$struc_method == "RFantibody" ~ "T",
      TRUE ~ input$ag_chain
    ) %>% unique()
    show_hchain <- case_when(
      input$struc_method %in% c("AbDockgen","RFantibody") ~ "H",
      TRUE ~ input$h_chain
    ) %>% unique()
    NGLVieweR(pdb_path) %>% 
      stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
      addRepresentation("cartoon",
                        param = list(sele = paste0(":",show_hchain), color = "red",
                                     opacity = 1)
      ) %>% 
      addRepresentation("surface", param = list(sele = paste0(":",show_agchain),
                                                color = "blue", opacity = 1))
  })
  
  output$struc_res_download <- downloadHandler(
    filename = function() {
      paste0(input_struc_data$pro_name,"_res.zip")
    },
    content = function(file) {
      zip_file <- paste0(sc_path,"/test/res/",input_struc_data$pro_name,"_res.zip")
      file.copy(zip_file, file)
    }
  )
} 

#options(shiny.fullstacktrace = TRUE)  # 显示完整堆栈跟踪
shinyApp(ui, server)


