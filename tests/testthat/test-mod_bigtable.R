test_that("module ui works", {
  ui <- mod_bigtable_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true(all(c("id", "height") %in% names(formals(mod_bigtable_ui))))
})

test_that("no data renders a placeholder table", {
  testServer(mod_bigtable_server, args = list(r_data = shiny::reactive(NULL)), {
    expect_null(full_data())
    expect_null(filtered())
    expect_false(is.null(output$table))
    expect_equal(output$mode_note, "")
  })
})

test_that("a small table stays client-side", {
  df <- data.frame(a = 1:20, b = letters[1:20], stringsAsFactors = FALSE)
  testServer(mod_bigtable_server, args = list(r_data = shiny::reactive(df)), {
    session$setInputs(search = "", sort_by = "", sort_dir = "asc", page_size = "25")
    expect_false(server_mode())
    expect_equal(nrow(page_data()), 20)
    expect_match(output$mode_note, "Client-side mode")
  })
})

test_that("a large table switches to server-side paging", {
  df <- data.frame(a = 1:500, b = paste0("row_", 1:500), stringsAsFactors = FALSE)
  testServer(mod_bigtable_server,
             args = list(r_data = shiny::reactive(df), server_threshold = 100), {
    session$setInputs(search = "", sort_by = "", sort_dir = "asc", page_size = "25")
    expect_true(server_mode())
    # Only the current page reaches the browser.
    expect_equal(nrow(page_data()), 25)
    expect_equal(page_data()$a, 1:25)
    expect_equal(n_pages(), 20)
    expect_match(output$mode_note, "Server-side mode")

    session$setInputs(next_page = 1)
    expect_equal(page_data()$a, 26:50)
    session$setInputs(next_page = 2)
    expect_equal(page_data()$a, 51:75)
    session$setInputs(prev_page = 1)
    expect_equal(page_data()$a, 26:50)
    expect_false(is.null(output$pager))
  })
})

test_that("search and sort apply to the whole table, not just the page", {
  df <- data.frame(a = 1:500, b = paste0("row_", 1:500), stringsAsFactors = FALSE)
  testServer(mod_bigtable_server,
             args = list(r_data = shiny::reactive(df), server_threshold = 100), {
    session$setInputs(search = "", sort_by = "a", sort_dir = "desc",
                      page_size = "10")
    expect_equal(page_data()$a, 500:491)

    # A match that lives far beyond the first page must still be found.
    session$setInputs(search = "row_499", sort_by = "", sort_dir = "asc")
    expect_equal(nrow(filtered()), 1)
    expect_equal(filtered()$a, 499)

    session$setInputs(search = "nothing_matches_this")
    expect_equal(nrow(filtered()), 0)
    expect_false(is.null(output$table))
  })
})

test_that("search treats special characters as literal text", {
  df <- data.frame(value = c("[", "plain"), stringsAsFactors = FALSE)
  testServer(mod_bigtable_server, args = list(r_data = shiny::reactive(df)), {
    session$setInputs(search = "[", sort_by = "", sort_dir = "asc", page_size = "25")
    expect_equal(filtered()$value, "[")
  })
})

test_that("changing the search resets to the first page", {
  df <- data.frame(a = 1:500, b = paste0("row_", 1:500), stringsAsFactors = FALSE)
  testServer(mod_bigtable_server,
             args = list(r_data = shiny::reactive(df), server_threshold = 100), {
    session$setInputs(search = "", sort_by = "", sort_dir = "asc", page_size = "25")
    session$setInputs(next_page = 1)
    session$setInputs(next_page = 2)
    expect_equal(page_data()$a[[1]], 51)
    session$setInputs(search = "row_")
    expect_equal(page_data()$a[[1]], 1)
  })
})

test_that("column selection restricts both display and search", {
  df <- data.frame(a = 1:10, b = letters[1:10],
                   secret = paste0("hidden_token_", 1:10),
                   stringsAsFactors = FALSE)
  testServer(mod_bigtable_server, args = list(r_data = shiny::reactive(df)), {
    session$setInputs(search = "", sort_by = "", sort_dir = "asc",
                      page_size = "25", columns = c("a", "b"))
    expect_equal(colnames(page_data()), c("a", "b"))
    # A term that only appears in a hidden column must not match.
    session$setInputs(search = "hidden_token")
    expect_equal(nrow(filtered()), 0)
    session$setInputs(columns = c("a", "b", "secret"))
    expect_equal(nrow(filtered()), 10)
  })
})

test_that("downloads contain the full search result", {
  df <- data.frame(a = 1:500, b = paste0("row_", 1:500), stringsAsFactors = FALSE)
  testServer(mod_bigtable_server,
             args = list(r_data = shiny::reactive(df), server_threshold = 100,
                         filename = "things"), {
    session$setInputs(search = "", sort_by = "", sort_dir = "asc", page_size = "25")
    f <- tempfile(fileext = ".tsv")
    write_tsv(f)
    written <- utils::read.delim(f)
    # Not just the 25 visible rows.
    expect_equal(nrow(written), 500)
    expect_match(name_tsv(), "^things_")
    unlink(f)

    csv <- tempfile(fileext = ".csv")
    write_csv(csv)
    expect_equal(nrow(utils::read.csv(csv)), 500)
    unlink(csv)
  })
})

test_that("a reactive filename is resolved at download time", {
  df <- data.frame(a = 1:3)
  testServer(mod_bigtable_server,
             args = list(r_data = shiny::reactive(df),
                         filename = shiny::reactive("dynamic")), {
    expect_match(name_csv(), "^dynamic_")
  })
})

test_that("a zero-row input does not error", {
  df <- data.frame(a = numeric(0), b = character(0), stringsAsFactors = FALSE)
  testServer(mod_bigtable_server, args = list(r_data = shiny::reactive(df)), {
    session$setInputs(search = "", sort_by = "", sort_dir = "asc", page_size = "25")
    expect_equal(nrow(page_data()), 0)
    expect_false(is.null(output$table))
    expect_equal(n_pages(), 1)
  })
})
