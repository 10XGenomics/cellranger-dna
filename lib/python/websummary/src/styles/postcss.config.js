// Copied from Bootstrap postcss config and prettiered
// https://github.com/twbs/bootstrap/blob/v4.0.0/build/postcss.config.js
'use strict';

module.exports = ctx => ({
  map: ctx.file.dirname.includes('examples')
    ? false
    : {
        inline: false,
        annotation: true,
        sourcesContent: true,
      },
  plugins: {
    autoprefixer: {
      cascade: false,
    },
  },
});
