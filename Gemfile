source "https://rubygems.org"

gem "github-pages", group: :jekyll_plugins

# Local copy of the theme, used only by _config_dev.yml so that local builds
# skip the 26 MB remote_theme download that jekyll-remote-theme repeats on
# every build. GitHub Pages ignores this Gemfile and keeps using remote_theme.
gem "minimal-mistakes-jekyll"

gem "liquid", "~> 4.0.4"
gem "csv"
gem "bigdecimal"
gem "webrick"
gem "tzinfo-data"
gem "wdm", "~> 0.2.0" if Gem.win_platform?

# If you have any plugins, put them here!
group :jekyll_plugins do
  gem "jekyll-paginate"
  gem "jekyll-sitemap"
  gem "jekyll-gist"
  gem "jekyll-feed"
  gem "jemoji"
  gem "jekyll-include-cache"
  gem "jekyll-algolia"
end
