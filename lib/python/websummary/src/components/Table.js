import React from 'react';
import map from 'lodash/map';

const Table = ({ header, rows }) => {
  return (
    <table className="table table-hover">
      {header && (
        <thead>
          <tr>{map(header, cell => <th>{cell}</th>)}</tr>
        </thead>
      )}
      <tbody>
        {map(rows, (row, i) => (
          <tr key={i}>{map(row, cell => <td key={cell}>{cell}</td>)}</tr>
        ))}
      </tbody>
    </table>
  );
};

export default Table;
