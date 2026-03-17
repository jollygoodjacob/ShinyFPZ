library(shiny)
library(sf)
library(dplyr)
library(cluster)
library(ggplot2)
library(FactoMineR)
library(RColorBrewer)
library(DT)
library(tidyr)
library(leaflet)
library(geojsonsf)

ui <- fluidPage(
  titlePanel("ShinyFPZ"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1) Upload data"),
      helpText("Upload a zipped shapefile containing .shp, .shx, .dbf, and .prj files."),
      fileInput("zipfile", "Upload zipped shapefile (.zip)", accept = ".zip"),
      
      hr(),
      h4("2) Clustering settings"),
      uiOutput("var_select_ui"),
      
      selectInput(
        "geo_var",
        "Categorical geology variable",
        choices = NULL,
        selected = "GEO"
      ),
      
      checkboxInput("scale_vars", "Min-max scale numeric variables", value = TRUE),
      checkboxInput("drop_na", "Drop rows with missing values", value = TRUE),
      
      radioButtons(
        "cluster_method",
        "Choose FPZ method",
        choices = c(
          "Manual k from silhouette curve" = "silhouette_manual",
          "Height cutoff" = "cutoff"
        ),
        selected = "silhouette_manual"
      ),
      
      conditionalPanel(
        condition = "input.cluster_method == 'silhouette_manual'",
        tagList(
          sliderInput(
            "k_range",
            "Range of k for silhouette curve",
            min = 2, max = 20, value = c(2, 10), step = 1
          ),
          uiOutput("k_pick_ui")
        )
      ),
      
      conditionalPanel(
        condition = "input.cluster_method == 'cutoff'",
        sliderInput(
          "cut_prop",
          "Cut tree at proportion of max height",
          min = 0.01, max = 1, value = 0.2, step = 0.01
        )
      ),
      
      actionButton("run", "Run FPZ analysis"),
      
      hr(),
      h4("3) Export outputs"),
      tagList(
        downloadButton("download_sil_plot", "Silhouette plot (.png)"),
        tags$br(), tags$br(),
        downloadButton("download_dend_plot", "Dendrogram (.png)"),
        tags$br(), tags$br(),
        downloadButton("download_pca_plot", "PCA biplot (.png)"),
        tags$br(), tags$br(),
        downloadButton("download_sil_table", "Silhouette table (.csv)"),
        tags$br(), tags$br(),
        downloadButton("download_cluster_table", "Cluster summary (.csv)"),
        tags$br(), tags$br(),
        downloadButton("download_gpkg", "Spatial output (.gpkg)")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("Preview", DTOutput("preview_table")),
        
        tabPanel(
          "Silhouette",
          verbatimTextOutput("best_k_text"),
          plotOutput("sil_plot", height = "500px"),
          br(),
          DTOutput("sil_table")
        ),
        
        tabPanel("Dendrogram", plotOutput("dend_plot", height = "500px")),
        tabPanel("PCA Biplot", plotOutput("pca_plot", height = "650px")),
        tabPanel("Cluster Summary", DTOutput("cluster_table")),
        
        tabPanel(
          "Leaflet Map",
          uiOutput("map_status"),
          leafletOutput("fpz_map", height = "700px"),
          br(),
          fluidRow(
            column(
              width = 6,
              div(
                style = "border: 1px solid #ddd; border-radius: 6px; padding: 10px; min-height: 360px;",
                h4("Map Diagnostics"),
                verbatimTextOutput("map_diag")
              )
            ),
            column(
              width = 6,
              div(
                style = "border: 1px solid #ddd; border-radius: 6px; padding: 10px; min-height: 360px;",
                h4("PCA Biplot"),
                plotOutput("pca_plot_map", height = "300px")
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  shp_data <- reactive({
    req(input$zipfile)
    
    td <- tempdir()
    unzip_dir <- file.path(td, paste0("shp_", as.integer(Sys.time())))
    dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
    unzip(input$zipfile$datapath, exdir = unzip_dir)
    
    shp_files <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)
    
    validate(
      need(length(shp_files) > 0, "No .shp file found inside the zip.")
    )
    
    sf::st_read(shp_files[1], quiet = TRUE)
  })
  
  observeEvent(shp_data(), {
    dat <- shp_data()
    nms <- names(dat)
    
    default_vars <- intersect(
      c("ELE", "PRE", "GEO", "VW", "VFW", "RAT", "RVS", "LVS", "DVS", "SIN"),
      nms
    )
    
    output$var_select_ui <- renderUI({
      selectInput(
        "vars",
        "Variables for FPZ clustering",
        choices = nms,
        selected = default_vars,
        multiple = TRUE
      )
    })
    
    updateSelectInput(
      session,
      "geo_var",
      choices = nms,
      selected = if ("GEO" %in% nms) "GEO" else nms[1]
    )
  })
  
  output$preview_table <- renderDT({
    req(shp_data())
    datatable(
      sf::st_drop_geometry(shp_data()),
      options = list(scrollX = TRUE, pageLength = 8),
      rownames = FALSE
    )
  })
  
  base_clustering <- eventReactive(input$run, {
    req(shp_data(), input$vars)
    
    center_points <- shp_data()
    
    validate(
      need(length(input$vars) >= 2, "Please choose at least 2 variables.")
    )
    
    dat0 <- center_points %>%
      dplyr::mutate(.row_id = dplyr::row_number()) %>%
      dplyr::select(.row_id, dplyr::all_of(input$vars))
    
    dat_tbl <- sf::st_drop_geometry(dat0)
    
    validate(
      need(input$geo_var %in% names(dat_tbl),
           "Selected geology variable is not in the chosen variables.")
    )
    
    dat_work <- dat_tbl
    dat_work[[input$geo_var]] <- as.factor(dat_work[[input$geo_var]])
    
    if (isTRUE(input$drop_na)) {
      dat_work <- tidyr::drop_na(dat_work)
    }
    
    validate(
      need(nrow(dat_work) >= 3, "Not enough complete rows remaining after filtering.")
    )
    
    vars_clust <- dat_work
    
    numeric_cols <- names(vars_clust)[sapply(vars_clust, is.numeric)]
    numeric_cols <- setdiff(numeric_cols, ".row_id")
    
    if (isTRUE(input$scale_vars) && length(numeric_cols) > 0) {
      vars_clust <- vars_clust %>%
        dplyr::mutate(
          dplyr::across(
            dplyr::all_of(numeric_cols),
            ~ {
              rng <- max(., na.rm = TRUE) - min(., na.rm = TRUE)
              if (is.na(rng) || rng == 0) 0 else (. - min(., na.rm = TRUE)) / rng
            }
          )
        )
    }
    
    diss <- cluster::daisy(
      vars_clust %>% dplyr::select(-.row_id),
      metric = "gower"
    )
    
    hc <- hclust(diss, method = "ward.D2")
    
    sil_df <- NULL
    suggested_k <- NULL
    
    if (input$cluster_method == "silhouette_manual") {
      ks <- seq(input$k_range[1], input$k_range[2])
      
      sil_avg <- sapply(ks, function(k) {
        cl <- cutree(hc, k = k)
        mean(cluster::silhouette(cl, diss)[, 3])
      })
      
      sil_df <- data.frame(
        k = ks,
        avg_sil = sil_avg
      )
      
      suggested_k <- ks[which.max(sil_avg)]
    }
    
    list(
      center_points = center_points,
      vars_clust = vars_clust,
      diss = diss,
      hc = hc,
      sil_df = sil_df,
      suggested_k = suggested_k
    )
  })
  
  output$k_pick_ui <- renderUI({
    req(base_clustering())
    req(input$cluster_method == "silhouette_manual")
    req(base_clustering()$sil_df)
    
    ks <- base_clustering()$sil_df$k
    suggested <- base_clustering()$suggested_k
    
    numericInput(
      "k_manual",
      "Choose number of FPZs (based on silhouette inflection/plateau)",
      value = suggested,
      min = min(ks),
      max = max(ks),
      step = 1
    )
  })
  
  analysis <- reactive({
    req(base_clustering())
    
    base <- base_clustering()
    center_points <- base$center_points
    vars_clust <- base$vars_clust
    diss <- base$diss
    hc <- base$hc
    
    if (input$cluster_method == "silhouette_manual") {
      req(input$k_manual)
      
      validate(
        need(
          input$k_manual >= min(base$sil_df$k) &&
            input$k_manual <= max(base$sil_df$k),
          "Selected k is outside the silhouette search range."
        )
      )
      
      FPZs <- cutree(hc, k = input$k_manual)
      sil <- cluster::silhouette(FPZs, diss)
      k_best <- input$k_manual
      sil_df <- base$sil_df
      
    } else {
      cutoff <- input$cut_prop * max(hc$height)
      FPZs <- cutree(hc, h = cutoff)
      sil <- cluster::silhouette(FPZs, diss)
      k_best <- length(unique(FPZs))
      sil_df <- NULL
    }
    
    vars_clust$FPZ <- factor(FPZs)
    
    vars_pca <- vars_clust %>%
      dplyr::mutate(
        !!input$geo_var := as.numeric(as.factor(.data[[input$geo_var]]))
      )
    
    pca_input <- vars_pca %>% dplyr::select(-.row_id, -FPZ)
    
    validate(
      need(ncol(pca_input) >= 2, "Need at least 2 variables for PCA.")
    )
    
    pca_result <- FactoMineR::PCA(
      pca_input,
      scale.unit = TRUE,
      graph = FALSE
    )
    
    percent_variance <- pca_result$eig[, 2]
    pc1_var <- round(percent_variance[1], 1)
    pc2_var <- round(percent_variance[2], 1)
    
    pca_coords <- as.data.frame(pca_result$ind$coord)
    pca_coords$FPZ <- vars_pca$FPZ
    
    var_coords <- as.data.frame(pca_result$var$coord)
    var_coords$varname <- rownames(var_coords)
    
    arrow_scale <- min(
      (max(pca_coords$Dim.1) - min(pca_coords$Dim.1)) /
        (max(var_coords$Dim.1) - min(var_coords$Dim.1)),
      (max(pca_coords$Dim.2) - min(pca_coords$Dim.2)) /
        (max(var_coords$Dim.2) - min(var_coords$Dim.2))
    ) * 0.8
    
    var_coords_scaled <- var_coords
    var_coords_scaled$Dim.1 <- var_coords$Dim.1 * arrow_scale
    var_coords_scaled$Dim.2 <- var_coords$Dim.2 * arrow_scale
    
    center_points_out <- center_points %>%
      dplyr::mutate(.row_id = dplyr::row_number()) %>%
      dplyr::left_join(
        vars_clust %>% dplyr::select(.row_id, FPZ),
        by = ".row_id"
      )
    
    cluster_summary <- vars_clust %>%
      dplyr::count(FPZ) %>%
      dplyr::rename(n_segments = n)
    
    list(
      center_points_out = center_points_out,
      hc = hc,
      diss = diss,
      sil = sil,
      sil_df = sil_df,
      k_best = k_best,
      pca_coords = pca_coords,
      var_coords_scaled = var_coords_scaled,
      pc1_var = pc1_var,
      pc2_var = pc2_var,
      cluster_summary = cluster_summary
    )
  })
  
  map_data <- reactive({
    req(analysis())
    
    shp <- analysis()$center_points_out
    
    validate(
      need(inherits(shp, "sf"), "Mapped output is not an sf object."),
      need(nrow(shp) > 0, "No features available to map."),
      need("FPZ" %in% names(shp), "Run FPZ analysis first to generate FPZ classes."),
      need(!all(sf::st_is_empty(shp)), "All geometries are empty."),
      need(!is.na(sf::st_crs(shp)), "Input layer has no CRS.")
    )
    
    shp <- shp[!sf::st_is_empty(shp), ]
    
    validate(
      need(nrow(shp) > 0, "No non-empty geometries available to map.")
    )
    
    shp <- sf::st_transform(shp, 4326)
    
    fpz_lab <- as.character(shp$FPZ)
    fpz_lab[is.na(fpz_lab)] <- "Unclassified"
    shp$FPZ_lab <- fpz_lab
    
    shp
  })
  
  sil_plot_obj <- reactive({
    req(analysis())
    
    if (!is.null(analysis()$sil_df)) {
      ggplot(analysis()$sil_df, aes(x = k, y = avg_sil)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_vline(xintercept = analysis()$k_best, linetype = "dashed") +
        geom_point(
          data = subset(analysis()$sil_df, k == analysis()$k_best),
          aes(x = k, y = avg_sil),
          size = 4
        ) +
        theme_bw() +
        labs(
          title = "Average silhouette width by number of FPZs",
          subtitle = "Choose k based on an elbow, inflection point, or stable plateau",
          x = "Number of FPZs (k)",
          y = "Average silhouette width"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title = element_text(face = "bold")
        )
    } else {
      NULL
    }
  })
  
  dend_plot_obj <- reactive({
    req(analysis())
    
    hc <- analysis()$hc
    k_best <- analysis()$k_best
    
    function() {
      op <- par(no.readonly = TRUE)
      on.exit(par(op))
      
      par(
        mar = c(4, 4, 3, 2) + 0.1,
        lwd = 1.3,
        cex.main = 1.1,
        cex.lab = 1,
        cex.axis = 0.9
      )
      
      plot(
        hc,
        main = "Hierarchical clustering dendrogram",
        xlab = "",
        sub = "",
        ylab = "Height",
        hang = -1,
        labels = FALSE
      )
      
      if (input$cluster_method == "cutoff") {
        cutoff_height <- input$cut_prop * max(hc$height)
        
        abline(
          h = cutoff_height,
          col = "red",
          lty = 2,
          lwd = 2
        )
        
        usr <- par("usr")
        text(
          x = usr[1] + 0.02 * diff(usr[1:2]),
          y = cutoff_height,
          labels = paste0("Cutoff = ", round(cutoff_height, 3)),
          pos = 3,
          col = "red",
          cex = 0.9
        )
      }
      
      rect.hclust(
        hc,
        k = k_best,
        border = "blue"
      )
    }
  })
  
  pca_plot_obj <- reactive({
    req(analysis())
    
    n_fpz <- length(unique(analysis()$pca_coords$FPZ))
    pal <- if (n_fpz <= 8) {
      RColorBrewer::brewer.pal(max(3, n_fpz), "Set1")[seq_len(n_fpz)]
    } else {
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_fpz)
    }
    
    ggplot(analysis()$pca_coords, aes(x = Dim.1, y = Dim.2, color = FPZ)) +
      geom_point(size = 2.5, shape = 19) +
      geom_segment(
        data = analysis()$var_coords_scaled,
        aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
        inherit.aes = FALSE,
        color = "black",
        arrow = arrow(length = grid::unit(0.2, "cm"))
      ) +
      geom_text(
        data = analysis()$var_coords_scaled,
        aes(x = Dim.1, y = Dim.2, label = varname),
        inherit.aes = FALSE,
        color = "black",
        vjust = -0.8,
        hjust = 0.5,
        size = 4
      ) +
      scale_color_manual(values = pal) +
      labs(
        title = "PCA - FPZ clustering on hydrogeomorphic features",
        x = paste0("PC1 (", analysis()$pc1_var, "% variance)"),
        y = paste0("PC2 (", analysis()$pc2_var, "% variance)"),
        color = "FPZ"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      )
  })
  
  output$best_k_text <- renderText({
    req(analysis())
    
    if (input$cluster_method == "silhouette_manual") {
      paste("User-selected number of FPZs:", analysis()$k_best)
    } else {
      paste("Number of FPZs from cutoff:", analysis()$k_best)
    }
  })
  
  output$sil_plot <- renderPlot({
    req(analysis())
    
    if (!is.null(analysis()$sil_df)) {
      print(sil_plot_obj())
    } else {
      plot(analysis()$sil, border = NA, main = "Silhouette plot")
    }
  })
  
  output$dend_plot <- renderPlot({
    req(dend_plot_obj())
    dend_plot_obj()()
  })
  
  output$pca_plot <- renderPlot({
    req(pca_plot_obj())
    print(pca_plot_obj())
  })
  
  output$pca_plot_map <- renderPlot({
    req(pca_plot_obj())
    print(
      pca_plot_obj() +
        guides(color = "none") +
        theme(
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8)
        )
    )
  })
  
  output$sil_table <- renderDT({
    req(analysis())
    
    if (!is.null(analysis()$sil_df)) {
      datatable(
        analysis()$sil_df,
        options = list(pageLength = 10, dom = "tip"),
        rownames = FALSE
      )
    }
  })
  
  output$cluster_table <- renderDT({
    req(analysis())
    datatable(
      analysis()$cluster_summary,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  output$map_status <- renderUI({
    if (is.null(input$run) || input$run < 1) {
      return(div(
        style = "margin-bottom:10px; color:#555;",
        "Run FPZ analysis to generate and map FPZ classes."
      ))
    }
    
    shp <- tryCatch(map_data(), error = function(e) NULL)
    
    if (is.null(shp)) {
      return(div(
        style = "margin-bottom:10px; color:#555;",
        "Map will appear after FPZ classification is successfully generated."
      ))
    }
    
    geom_type <- unique(as.character(sf::st_geometry_type(shp, by_geometry = TRUE)))
    n_cls <- length(unique(shp$FPZ_lab))
    
    div(
      style = "margin-bottom:10px; color:#555;",
      paste(
        "Mapped", nrow(shp), "features across", n_cls, "FPZ classes.",
        "Geometry:", paste(geom_type, collapse = ", ")
      )
    )
  })
  
  output$fpz_map <- renderLeaflet({
    shp <- map_data()
    
    geom_type <- unique(as.character(sf::st_geometry_type(shp, by_geometry = TRUE)))
    geom_type <- stats::na.omit(geom_type)
    
    is_point <- length(geom_type) > 0 && any(grepl("POINT", geom_type))
    is_line  <- length(geom_type) > 0 && any(grepl("LINESTRING", geom_type))
    
    fpz_vals <- sort(unique(as.character(shp$FPZ_lab)))
    
    base_cols <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
      "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#1B9E77"
    )
    
    cols <- rep(base_cols, length.out = length(fpz_vals))
    names(cols) <- fpz_vals
    
    if (is_point) {
      coords <- sf::st_coordinates(shp)
      
      map_df <- data.frame(
        lon = coords[, 1],
        lat = coords[, 2],
        FPZ = as.character(shp$FPZ_lab),
        stringsAsFactors = FALSE
      )
      
      map_df$col <- unname(cols[map_df$FPZ])
      
      center_lon <- mean(range(map_df$lon, na.rm = TRUE))
      center_lat <- mean(range(map_df$lat, na.rm = TRUE))
      
      leaflet(map_df) %>%
        addProviderTiles(providers$Esri.WorldTopoMap) %>%
        setView(lng = center_lon, lat = center_lat, zoom = 9) %>%
        addCircleMarkers(
          lng = ~lon,
          lat = ~lat,
          radius = 5,
          stroke = TRUE,
          color = ~col,
          weight = 2,
          opacity = 1,
          fill = FALSE,
          popup = ~paste("FPZ:", FPZ)
        ) %>%
        addLegend(
          position = "bottomright",
          colors = unname(cols),
          labels = unname(fpz_vals),
          title = "FPZ",
          opacity = 1
        )
      
    } else if (is_line) {
      bb <- sf::st_bbox(shp)
      center_lon <- mean(c(bb["xmin"], bb["xmax"]))
      center_lat <- mean(c(bb["ymin"], bb["ymax"]))
      
      fpz_vals <- sort(unique(as.character(shp$FPZ_lab)))
      
      base_cols <- c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62"
      )
      
      cols <- rep(base_cols, length.out = length(fpz_vals))
      names(cols) <- fpz_vals
      
      m <- leaflet() %>%
        addProviderTiles(providers$Esri.WorldImagery) %>%
        setView(lng = center_lon, lat = center_lat, zoom = 9)
      
      for (cls in fpz_vals) {
        shp_sub <- shp[shp$FPZ_lab == cls, ]
        shp_json <- geojsonsf::sf_geojson(shp_sub)
        this_col <- as.character(unname(cols[cls]))
        
        m <- m %>%
          addGeoJSON(
            shp_json,
            color = this_col,
            weight = 4,
            opacity = 0.9
          )
      }
      
      m
      
    } else {
      validate(
        need(FALSE, paste("Unsupported geometry type:", paste(geom_type, collapse = ", ")))
      )
    }
  })
  
  output$map_diag <- renderPrint({
    shp <- tryCatch(map_data(), error = function(e) NULL)
    if (is.null(shp)) return("Map not ready.")
    
    geom_type <- unique(as.character(sf::st_geometry_type(shp, by_geometry = TRUE)))
    
    list(
      n_features = nrow(shp),
      fpz_counts = table(shp$FPZ_lab, useNA = "ifany"),
      geometry_types = geom_type,
      any_na_geometry_type = any(is.na(geom_type)),
      bbox = sf::st_bbox(shp)
    )
  })
  
  output$download_sil_plot <- downloadHandler(
    filename = function() "silhouette_plot.png",
    content = function(file) {
      req(analysis())
      
      if (!is.null(analysis()$sil_df)) {
        ggplot2::ggsave(file, plot = sil_plot_obj(), width = 8, height = 5, dpi = 300)
      } else {
        png(file, width = 8, height = 5, units = "in", res = 300)
        plot(analysis()$sil, border = NA, main = "Silhouette plot")
        dev.off()
      }
    }
  )
  
  output$download_dend_plot <- downloadHandler(
    filename = function() "dendrogram.png",
    content = function(file) {
      req(dend_plot_obj())
      png(file, width = 8, height = 5, units = "in", res = 300)
      dend_plot_obj()()
      dev.off()
    }
  )
  
  output$download_pca_plot <- downloadHandler(
    filename = function() "pca_biplot.png",
    content = function(file) {
      req(pca_plot_obj())
      ggplot2::ggsave(file, plot = pca_plot_obj(), width = 8, height = 6.5, dpi = 300)
    }
  )
  
  output$download_sil_table <- downloadHandler(
    filename = function() "silhouette_table.csv",
    content = function(file) {
      req(analysis())
      
      if (!is.null(analysis()$sil_df)) {
        write.csv(analysis()$sil_df, file, row.names = FALSE)
      } else {
        write.csv(
          data.frame(message = "Silhouette table not available for cutoff mode"),
          file,
          row.names = FALSE
        )
      }
    }
  )
  
  output$download_cluster_table <- downloadHandler(
    filename = function() "cluster_summary.csv",
    content = function(file) {
      req(analysis())
      write.csv(analysis()$cluster_summary, file, row.names = FALSE)
    }
  )
  
  output$download_gpkg <- downloadHandler(
    filename = function() "FPZ_output.gpkg",
    content = function(file) {
      req(analysis())
      sf::st_write(analysis()$center_points_out, file, delete_dsn = TRUE, quiet = TRUE)
    }
  )
}

shinyApp(ui, server)