import React, { Component } from 'react';
import Plotly from './Plotly';
import map from 'lodash/map';
import merge from 'lodash/merge';
import times from 'lodash/times';

export default class Plot extends Component {
  componentDidMount() {
    const { data, layout = {}, config = {} } = this.props;

    const defaultLayout = {
      xaxis: { fixedrange: true },
      yaxis: { fixedrange: true },
      font: { family: 'DIN Next LT Pro' },
      // hovermode: 'closest',
      margin: { t: 24 },
    };

    const defaultConfig = {
      displaylogo: false,
      displayModeBar: false,
      showLink: false,
      staticPlot: true,
    };

    /*
    const mergedData = map(data, datum =>
      merge(datum, {
        hoverinfo: 'text',
        text: map(times(datum.y.length),
          i =>
            `${datum.x[i].toLocaleString('en-US')} ${
              layout.xaxis.title
            }<br>${datum.y[i].toLocaleString('en-US')} ${layout.yaxis.title}`
        ),
      })
    );
    */
    const mergedLayout = merge({}, defaultLayout, layout);
    const mergedConfig = merge({}, defaultConfig, config);

    Plotly.react(this.plot, {
      data,
      layout: mergedLayout,
      config: mergedConfig,
    });

    this.resizeHandler = () => Plotly.Plots.resize(this.plot);
    window.addEventListener('resize', this.resizeHandler);
  }

  componentWillUnmount() {
    window.removeEventListener('resize', this.resizeHandler);
  }

  render() {
    return (
      <div
        className="plot-wrapper"
        ref={plot => {
          this.plot = plot;
        }}
      />
    );
  }
}
