import Plotly from 'plotly.js/lib/core';

Plotly.register([
    require('plotly.js/lib/bar'),
    require('plotly.js/lib/scatter'),
    require('plotly.js/lib/sankey'),
]);

export default Plotly;
