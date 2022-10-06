mod_jb2_ui <- function(id) {
  tagList(
    # this adds to the browser to the UI, and specifies the output ID in the server
    JBrowseROutput(NS(id,"browserOutput"))
  )
}

mod_jb2_server <- function(id) {
  data_server <- serve_data("./JBrowser_data/", port = 5000)
  
  moduleServer(id, function(input, output, session) {
    # create the necessary JB2 assembly configuration
    assembly <- assembly(
      "http://127.0.0.1:5000/genome.fasta.gz",
      bgzip = TRUE
    )
    # create configuration for a JB2 GFF Feature Track
    annotations_track <- track_feature(
      # "https://jbrowse.org/genomes/sars-cov2/sars-cov2-annotations.sorted.gff.gz",
      "http://127.0.0.1:5000/trnascan_out.gff.gz",
      assembly
    )
    annotations_track2 <- track_feature(
      # "https://jbrowse.org/genomes/sars-cov2/sars-cov2-annotations.sorted.gff.gz",
      "http://127.0.0.1:5000/vpf.genbank.sorted.gff.gz",
      assembly
    )
    # create the tracks array to pass to browser
    tracks <- tracks(annotations_track, annotations_track2)
    theme <- theme("#333", "#ff6200")
    default_session <- default_session(
      assembly,
      tracks(annotations_track2)
    )
    # link the UI with the browser widget
    output$browserOutput <- renderJBrowseR(
      JBrowseR(
        "View",
        assembly = assembly,
        tracks = tracks,
        theme = theme,
        defaultSession = default_session
      )
    )
    session$onSessionEnded(function() {
      data_server$stop_server()
    })
  })
}
