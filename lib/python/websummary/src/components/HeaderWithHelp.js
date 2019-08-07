import React, { PureComponent, Fragment } from 'react';
import QuestionCircle from './icons/QuestionCircle';

export default class TitleWithHelp extends PureComponent {
  state = { helpShown: false };

  toggleHelp() {
    const { helpShown } = this.state;
    this.setState({ helpShown: !helpShown });
  }

  render() {
    const { title, helpText } = this.props;
    const { helpShown } = this.state;

    const helpDisplay = helpShown ? 'block' : 'none';

    return (
      <Fragment>
        <div className="header-wrapper">
          <div>{title && <h2 id={title}>{title}</h2>}</div>
          <div>
            {helpText && (
              <QuestionCircle
                onClick={() => this.toggleHelp()}
                className="header-questionMark"
              />
            )}
          </div>
        </div>
        <div
          style={{ display: helpDisplay }}
          dangerouslySetInnerHTML={{ __html: helpText }}
        />
      </Fragment>
    );
  }
}
