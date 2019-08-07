const path = require('path');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
const UglifyJsPlugin = require('uglifyjs-webpack-plugin');
const OptimizeCSSAssetsPlugin = require('optimize-css-assets-webpack-plugin');

module.exports = {
  entry: './index.js',
  output: {
    path: path.resolve(__dirname, '../dist'),
    filename: 'tenx-websummary-script.min.js',
  },
  module: {
    rules: [
      { test: /\.js$/, loader: 'babel-loader', exclude: /node_modules/ },
      // for plotly https://github.com/plotly/plotly-webpack#explanations
      { test: /\.js$/, use: ['ify-loader'] },
      // Mostly copied from the docs; removed style-loader because we dont want
      // CSS in the page and instead use ExtractTextPlugin to generate a file
      // https://getbootstrap.com/docs/4.0/getting-started/webpack/#importing-precompiled-sass
      {
        test: /\.(scss)$/,
        use: [
          { loader: MiniCssExtractPlugin.loader },
          // translates CSS into CommonJS modules
          { loader: 'css-loader' },
          {
            loader: 'postcss-loader', // Run post css actions
            options: {
              // post css plugins, can be exported to postcss.config.js
              plugins: () => [require('precss'), require('autoprefixer')],
            },
          },
          // compiles Sass to CSS
          { loader: 'sass-loader' },
        ],
      },
      {
        test: /\.(jpe?g|png|ttf|eot|svg|woff(2)?|ico)$/,
        use: 'url-loader',
      },
    ],
  },
  plugins: [
    new MiniCssExtractPlugin({ filename: 'tenx-websummary-styles.min.css' }),
  ],
  optimization: {
    minimizer: [
      new UglifyJsPlugin({ cache: true, parallel: true }),
      new OptimizeCSSAssetsPlugin({}),
    ],
  },
};
