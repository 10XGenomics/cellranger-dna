import React, { PureComponent, Fragment } from 'react';
import TimesOctagon from './icons/TimesOctagon';
import ExclamationTriangle from './icons/ExclamationTriangle';
import map from 'lodash/map';

export default class Alert extends PureComponent {
  state = { alertTableShown: false };

  toggleAlertTable() {
    const { alertTableShown } = this.state;
    this.setState({ alertTableShown: !alertTableShown });
  }

  render() {
    const { alarms: alerts } = this.props;
    const { alertTableShown } = this.state;

    if (alerts.length === 0) {
      return null;
    }

    const warnings = alerts.filter(alert => alert.level === 'WARN');
    const errors = alerts.filter(alert => alert.level === 'ERROR');

    const problemElements = [];

    if (errors.length > 0) {
      let errorString = '';

      if (errors.length > 1) {
        errorString = `${errors.length} errors`;
      } else {
        errorString = `1 error`;
      }

      problemElements.push(<TimesOctagon key="danger-icon" />);
      problemElements.push(
        <span key="danger" className="text-brand-red">
          {errorString}
        </span>
      );
    }

    if (warnings.length > 0 && errors.length > 0) {
      problemElements.push(' and ');
    }

    if (warnings.length > 0) {
      let warningString = '';

      if (warnings.length > 1) {
        warningString = `${warnings.length} warnings`;
      } else {
        warningString = `1 warning`;
      }

      problemElements.push(<ExclamationTriangle key="warning-icon" />);
      problemElements.push(
        <span key="warning" className="text-brand-yellow">
          {warningString}
        </span>
      );
    }

    const alertTableDisplay = alertTableShown ? 'block' : 'none';

    return (
      <div className="container alert-container">
        <div>
          <h2 className="alert-header">Alerts</h2>
          <div className="alert-analysis">
            The analysis detected {problemElements}.
          </div>
          <table className="table table-hover table-sm">
            <thead>
              <tr>
                <td />
                <td>Alert</td>
                <td>Value</td>
                <td>Detail</td>
              </tr>
            </thead>
            <tbody>
              {map(
                errors,
                ({ title, formatted_value: formattedValue, message }, i) => {
                  return (
                    <tr key={i}>
                      <td>
                        <TimesOctagon />
                      </td>
                      <td>{title}</td>
                      <td>{formattedValue}</td>
                      <td>{message}</td>
                    </tr>
                  );
                }
              )}
              {map(
                warnings,
                ({ title, formatted_value: formattedValue, message }, i) => {
                  return (
                    <tr key={i}>
                      <td>
                        <ExclamationTriangle />
                      </td>
                      <td>{title}</td>
                      <td>{formattedValue}</td>
                      <td>{message}</td>
                    </tr>
                  );
                }
              )}
            </tbody>
          </table>
        </div>
      </div>
    );
  }
}
