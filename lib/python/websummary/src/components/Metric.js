import React, { Component } from 'react';

const Metric = ({ name, metric, threshold }) => {
  return (
    <div className="metric-wrapper">
      <div className={`metric-number metric-threshold-${threshold}`}>
        {metric}
      </div>
      <div className="metric-name">{name}</div>
    </div>
  );
};

export default Metric;
