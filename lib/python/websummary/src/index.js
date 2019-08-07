import Alert from './components/Alert';
import HeaderWithHelp from './components/HeaderWithHelp.js';
import Metric from './components/Metric';
import Namescription from './components/Namescription';
import Navbar from './components/Navbar';
import Plot from './components/Plot';
import TableMetric from './components/TableMetric';
import Table from './components/Table';

import React from 'react';
import ReactDOM from 'react-dom';

import map from 'lodash/map';
import get from 'lodash/get';

import favicon from './assets/favicon.ico';

// https://reactjs.org/docs/jsx-in-depth.html#choosing-the-type-at-runtime
const components = {
  Alert,
  Metric,
  Navbar,
  Plot,
  TableMetric,
  Table,
};

// add favicon link tag
const faviconElement = document.createElement('link');
faviconElement.setAttribute('href', `data:image/x-icon;base64,${favicon}`);
document.head.appendChild(faviconElement);

// stores titles from header components for rendering in the navbar
class TitlesManager {
  titles = [];

  addTitle(title) {
    this.titles.push(title);
    renderNavbar();
  }
  getTitles() {
    return titles;
  }
}
const titlesManager = new TitlesManager();

const sampleId = get(data, 'sample.id');
const sampleDescription = get(data, 'sample.description');
const pipelineName = get(data, 'sample.pipeline');

let namescription = sampleId;
if (sampleDescription) {
  namescription += ` - ${sampleDescription}`;
}
document.title = namescription;

// TODO stop using data.sample after pipeline params gets fixed
const renderNavbar = () => {
  ReactDOM.render(
    <Navbar pipelineName={pipelineName} titles={titlesManager.titles} />,
    document.querySelector('.navbar-wrapper')
  );
};
renderNavbar();

ReactDOM.render(
  <Namescription namescription={namescription} />,
  document.querySelector('.namescription-wrapper')
);

const alarms = data.alarms || { alarms: [] };
ReactDOM.render(
  <Alert {...alarms} />,
  document.querySelector('.alert-wrapper')
);

map(document.querySelectorAll('[data-header]'), node => {
  const { title } = node.dataset;

  const helpDom = node.querySelector('.helpText');
  const helpText = helpDom ? helpDom.innerHTML : null;

  titlesManager.addTitle(title);
  ReactDOM.render(<HeaderWithHelp title={title} helpText={helpText} />, node);
});

map(document.querySelectorAll('[data-component]'), node => {
  const { key, component: componentType } = node.dataset;

  const Component = components[componentType];
  const props = { ...data[key] };

  ReactDOM.render(<Component {...props} />, node);
});

// these styles get extracted into a different css file
import './styles/styles.scss';
