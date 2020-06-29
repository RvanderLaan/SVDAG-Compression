
import { QuadTree, NULL_NODE, canvasSize } from './main.js';

class GraphRenderer {
  constructor(
    public tree: QuadTree,
    public canvas: HTMLCanvasElement,
  ) {
    canvas.width = canvasSize;
    canvas.height = canvasSize;
  }
  render() {
    const { width, height } = this.canvas;

    const { nLevels } = this.tree;

    const nRows = nLevels * 2 - 1; // 2 rows per level: one for nodes, next for lines (but not for leaves)

    const rowHeight = height / nRows;

    const ctx = this.canvas.getContext('2d');
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, width, height);

    // Prepare for drawing nodes
    ctx.strokeStyle = 'black';
    ctx.lineWidth = 4;

    const nodeRadius = Math.min(
      Math.min(rowHeight / 3, 64),
      0.5 * width / (this.tree.data[nLevels - 1].length + 1),
    );
    
    let globalNodeCount = 0;
    for (let level = 0; level < nLevels; level++) {
      const nodes = this.tree.data[level];

      const y = level * 2 * rowHeight + rowHeight / 2;
      const childY = (level + 1) * 2 * rowHeight + rowHeight / 2;

      // Draw nodes of this level
      for (let i = 0; i < nodes.length; i++) {
        // e.g. node 2 out of 5:  0.5w + (2 - 2.5) * w = 0.5w - 0.5 w = 0;
        const x = (i + 1) / (nodes.length + 1) * width;
        console.log('Drawing node', x, y, nodeRadius);

        
        if (level !== nLevels - 1) {
          ctx.beginPath();
          ctx.arc(x, y, nodeRadius, 0, 2 * Math.PI);
          ctx.stroke();

          ctx.fillStyle = 'black';
          ctx.font = `${nodeRadius * 1.5}px Arial`;
          const name = String.fromCharCode('A'.charCodeAt(0) + globalNodeCount);
          const size = ctx.measureText(name);
          ctx.fillText(name, x - size.width / 2, y + nodeRadius * 0.5);

          globalNodeCount++;

          // Draw lines to children
          const childLevelNodes = this.tree.data[level + 1];
          for (const c of nodes[i].children) {
            if (c !== NULL_NODE) {
              const childX = (c + 1) / (childLevelNodes.length + 1) * width;
              ctx.beginPath();
              ctx.moveTo(x, y + nodeRadius);
              ctx.lineTo(childX, childY - nodeRadius);
              ctx.stroke();
            }
          }
        } else {
          // Leaves are squares with 4 pixels
          const nr = nodeRadius * 0.9;
          const dnr = nr * 2;
          ctx.strokeRect(x - nr, y - nr, dnr, dnr);
          for (let c = 0; c < 4; c++) {
            if (nodes[i].childrenBitmask[c]) {
              ctx.fillRect(x + (c % 2 - 1) * nr, y + Math.floor(c / 2 - 1) * nr, nr, nr);
            }
          }
        }
      }
    }
  }
}

export default GraphRenderer;
