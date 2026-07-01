# Load required libraries -------------------------------------------------

# httr is used to send HTTP requests to the Wikimedia API
library(httr)

# jsonlite is used to convert the JSON response from the API into R objects
library(jsonlite)

# dplyr is used here mainly for bind_rows(), which combines many data frames
library(dplyr)


# Set working directory ---------------------------------------------------

# This tells R where to look for input files and where to save output files.
# Change this path if the location of your files changes.
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Data/Wikipedia data")


# Load language codes -----------------------------------------------------

# This file should contain a column called "Code", with language codes such as:
# en, fr, de, es, it, ja, etc.
#
# These codes are used to build Wikimedia project names such as:
# en.wikipedia.org
# fr.wikipedia.org
# de.wikipedia.org
language_list <- read.csv("ISO Codes.csv")


# Define function to retrieve Wikipedia pageviews -------------------------

# This function retrieves the total number of monthly pageviews for one
# Wikipedia article in one language.
#
# Arguments:
# article = Wikipedia page title, usually as it appears in the URL
#           Example: "Acinonyx_jubatus"
#
# lang    = language code
#           Example: "en" for English, "fr" for French
#
# start   = first date to query, in YYYYMMDD format
#           Important: Wikimedia pageview data only start from 2015-07-01
#
# end     = final date to query, in YYYYMMDD format
#           The API returns data up to this date
#
# The function returns a data frame with:
# language = the language code queried
# views    = total pageviews over the requested period
get_pageviews <- function(article,
                          lang = "en",
                          start = "20050101",
                          end = "20250101") {
  
  # Print the language currently being queried.
  # This is useful because, if the loop crashes or is slow,
  # you can see which language caused the problem.
  print(lang)
  
  
  # Base URL for the Wikimedia Pageviews API.
  # This part is constant for all requests.
  base_url <- "https://wikimedia.org/api/rest_v1/metrics/pageviews/per-article"
  
  
  # Build the Wikimedia project name.
  #
  # For English Wikipedia, the project must be:
  # en.wikipedia.org
  #
  # For French Wikipedia:
  # fr.wikipedia.org
  #
  # The older/incomplete version "en.wikipedia" will not work reliably.
  project <- paste0(lang, ".wikipedia.org")
  
  
  # Construct the full API request URL.
  #
  # The expected API structure is:
  #
  # /metrics/pageviews/per-article/
  #   {project}/
  #   {access}/
  #   {agent}/
  #   {article}/
  #   {granularity}/
  #   {start}/
  #   {end}
  #
  # Here:
  # all-access = includes desktop, mobile web, and mobile app
  # user       = excludes spiders/crawlers as much as Wikimedia can identify them
  # monthly    = returns one value per month
  #
  # URLencode(article, reserved = TRUE) makes the article title safe for a URL.
  # This is useful if page titles contain spaces, accents, parentheses, etc.
  url <- paste0(
    base_url, "/",
    project, "/",
    "all-access/user/",
    URLencode(article, reserved = TRUE),
    "/monthly/",
    start, "/",
    end
  )
  
  
  # Send the request to the Wikimedia API.
  #
  # The user_agent() line is important because Wikimedia requests clients
  # to identify themselves. Replace the email with your real contact email
  # if you want to follow Wikimedia API etiquette properly.
  res <- GET(
    url,
    user_agent("ZooCancerWikipediaViews/1.0 (antoine.dujon@example.com) httr")
  )
  
  
  # Check whether the request was successful.
  #
  # HTTP status code 200 means success.
  #
  # Common non-200 codes:
  # 404 = page not found, often because the article title does not exist
  #       in that language
  # 403 = request forbidden, often linked to API access rules/user-agent
  # 429 = too many requests
  #
  # If the request failed, return NA for that language rather than stopping
  # the whole loop.
  if (status_code(res) != 200) {
    
    # Print a useful diagnostic message.
    # This helps identify which languages failed and why.
    message("Failed: ", lang, " | status: ", status_code(res))
    
    # Return a one-row data frame with NA views.
    return(data.frame(language = lang, views = NA))
  }
  
  
  # Convert the API response from JSON text into an R object.
  #
  # content(res, as = "text") extracts the raw JSON response as text.
  # fromJSON() turns that text into a list/data frame structure.
  dat <- fromJSON(content(res, as = "text", encoding = "UTF-8"))
  
  
  # Extract monthly views and sum them.
  #
  # dat$items$views contains one value per month.
  #
  # na.rm = TRUE means that any missing values are ignored.
  total_views <- sum(dat$items$views, na.rm = TRUE)
  
  
  # Return the result as a simple one-row data frame.
  #
  # This format is convenient because later we can bind many languages
  # together using bind_rows().
  data.frame(language = lang, views = total_views)
}


# Prepare list of languages -----------------------------------------------

# Extract the language codes from the language_list data frame.
#
# This assumes your CSV has a column named exactly "Code".
languages <- language_list$Code


# Define target article ----------------------------------------------------

# This is the Wikipedia page title to query.
#
# For English Wikipedia, the cheetah species page is:
# https://en.wikipedia.org/wiki/Acinonyx_jubatus
#
# Therefore the article title is:
# Acinonyx_jubatus
#
# Important limitation:
# This exact title may not exist in other languages.
# For example, the French page may be "Guépard", not "Acinonyx_jubatus".
# Therefore some languages may return NA even if the species has a page.
article <- "Macaca_sylvanus"


# Query all languages ------------------------------------------------------

# lapply() applies the get_pageviews() function to every language code.
#
# For each language, it runs:
# get_pageviews(article = article, lang = lang)
#
# The result is a list of small one-row data frames.
#
# Sys.sleep(0.1) pauses for 0.1 seconds between requests.
# This is polite and reduces the risk of being rate-limited by the API.
#
# bind_rows() combines all the one-row data frames into one large data frame.
results <- lapply(languages, function(lang) {
  
  # Pause briefly between API calls
  Sys.sleep(0.1)
  
  # Retrieve pageviews for this article-language combination
  get_pageviews(article = article, lang = lang)
  
}) %>%
  bind_rows()


# Inspect results ----------------------------------------------------------

# Print the full results table.
#
# Each row corresponds to one language.
# views is the summed number of monthly pageviews over the selected period.
print(results)


# Sum total views across languages ----------------------------------------

# Sum the views column across all languages.
#
# na.rm = TRUE means languages with failed requests or missing pages
# are ignored rather than causing the total to become NA.
sum(results$views, na.rm = TRUE)