/* This is a two column table where the right column is:
 * in a monospaced font
 * right aligned
 * commatized if it's a number
 * It requires that each row only has a length of 2.
 *
 * If the key has the string "path" in it, it assumes its a filesystem path
 * and does magic truncation on it which un-truncates on click.
 */
import React, { PureComponent } from 'react';
import map from 'lodash/map';

export default class TableMetric extends PureComponent {
  render() {
    const { rows } = this.props;

    return (
      <table className="table table-hover table-sm">
        <tbody>
          {map(rows, (row, i) => {
            if (row.length !== 2) {
              console.log(`TableMetric: Row ${i} does not have two elements`);
            }

            const [key, val] = row;
            const isPath = key.toLowerCase().includes('path');

            return (
              <tr key={i}>
                <td>{key}</td>
                {!isPath && <td className="tableMetric-value">{val}</td>}
                {isPath && <TruncatedMetric text={val} />}
              </tr>
            );
          })}
        </tbody>
      </table>
    );
  }
}

class TruncatedMetric extends PureComponent {
  state = { truncated: true };

  render() {
    const { text } = this.props;
    const { truncated } = this.state;

    return (
      <td
        className="tableMetric-value"
        onMouseDown={() => this.setState({ truncated: false })}
      >
        {truncated ? (
          <div className="tableMetric-truncate">&lrm;{`${text}`}&lrm;</div>
        ) : (
          <div className="tableMetric-fixedWidth">{text}</div>
        )}
      </td>
    );
  }
}
