
const wait = async (time = 500) => new Promise((resolve) => setTimeout(resolve, time));

/////////////////////
///// QUAD TREE /////
/////////////////////

const nullNode = Number.MAX_VALUE;

const NCHILDREN = 4; // quad tree

type NodeData = QuadTreeNode[][];

class QuadTreeNode {
  childrenBitmask: boolean[];
  children: number[];
  outsideMask: boolean[];

  // DEBUG
  _x: number;
  _y: number;
  _size: number;
  _level: number;


  constructor(x: number, y: number, size: number, level: number) {
    this.childrenBitmask = []
    this.children = [];
    this.outsideMask = [];
    for (let i = 0; i < NCHILDREN; i++) {
      this.childrenBitmask[i] = false;
      this.children[i] = nullNode;
      this.outsideMask[i] = false;
    }

    this._x = x;
    this._y = y;
    this._size = size;
    this._level = level;
  }

  setChild(i: number, value: number) {
    this.childrenBitmask[i] = true;
    this.children[i] = value;
  }

  unsetChild(i: number) {
    this.childrenBitmask[i] = false;
    this.children[i] = nullNode;
  }

  hasChildren(): boolean {
    for (let i = 0; i < 4; i++) {
      if (this.childrenBitmask[i]) {
        return true;
      }
    }
    return false;
  }

  isInside(): boolean {
    for (let i = 0; i < 4; i++) {
      if (!this.outsideMask[i]) {
        return false;
      }
    }
    return true;
  }

  async draw(ctx: CanvasRenderingContext2D, x: number, y: number, size: number, lev: number, data: NodeData, nbs?: QuadTreeNode[]) {
    ctx.strokeStyle = 'green';
    ctx.strokeRect(x, y, size, size);
    // await wait(100);

    const hs = size / 2; // half size
    if (lev === data.length - 1) {
      ctx.fillStyle = 'black';
      for (let i = 0; i < NCHILDREN; i++) {
        if (this.childrenBitmask[i]) {
          ctx.fillRect(x + (i % 2) * hs, y + Math.floor(i / 2) * hs, hs, hs);
        }
      }
    } else {
      for (let i = 0; i < NCHILDREN; i++) {
        if (this.childrenBitmask[i]) {
          await data[lev + 1][this.children[i]].draw(ctx, x + (i % 2) * hs, y + Math.floor(i / 2) * hs, hs, lev + 1, data);
        }
      }
    }

    ctx.fillStyle = 'rgba(150, 150, 250, 0.2)';
    ctx.strokeStyle = 'rgba(150, 150, 250, 0.5)';
    ctx.lineWidth = 4;
    for (let i = 0; i < NCHILDREN; i++) {
      if (this.outsideMask[i]) {
        ctx.strokeRect(x + (i % 2) * hs, y + Math.floor(i / 2) * hs, hs, hs);
        ctx.fillRect(x + (i % 2) * hs, y + Math.floor(i / 2) * hs, hs, hs);
      }
    }
    ctx.lineWidth = 1;

    ctx.fillStyle = 'rgba(250, 150, 150, 0.5)';
    for (let i = 0; i < NCHILDREN; i++) {
      if (nbs && this.childrenBitmask[i] && nbs.some((nb) => nb === data[lev + 1][this.children[i]])) {
        ctx.fillRect(x + (i % 2) * hs, y + Math.floor(i / 2) * hs, hs, hs);
      }
    }
  }

  drawDirect(ctx: CanvasRenderingContext2D, lineWidth = 4, strokeStyle = 'rgba(150, 150, 250, 0.5)') {
    ctx.strokeStyle = strokeStyle;
    ctx.lineWidth = lineWidth;
    ctx.strokeRect(this._x, this._y, this._size, this._size);
  }
}

class QuadTree {
  nLevels: number;
  data: NodeData;

  constructor(nLevels: number) {
    this.nLevels = nLevels;
    this.data = [];
    for (let i = 0; i < nLevels; i++) {
      this.data[i] = [];
    }
    this.data[0].push(new QuadTreeNode(0, 0, canvasSize, 0)); // root node
  }

  setVoxel(x: number, y: number) {
    const r = Math.pow(2, this.nLevels);
    let parent = 0; // start with root
    let parentStartX = 0, parentStartY = 0;
    for (let lev = 0; lev < this.nLevels; lev++) {
      const n = this.data[lev][parent];
      const childCellSize = Math.pow(2, this.nLevels - lev - 1); // lev 0 -> 8, lev 1 -> 4, lev 2 -> 2, ...
      let childIndex = 0;
      if (x >= parentStartX + childCellSize) {
        childIndex += 1;
        parentStartX += childCellSize;
      }
      if (y >= parentStartY + childCellSize) {
        childIndex += 2;
        parentStartY += childCellSize;
      }
      if (lev < nLevels - 1) {
        if (!n.childrenBitmask[childIndex]) {
          const newParent = this.data[lev + 1].length;
          n.setChild(childIndex, newParent);
          const hs = n._size/2;
          this.data[lev + 1].push(new QuadTreeNode(
            n._x + hs * (childIndex % 2),
            n._y + hs * Math.floor(childIndex / 2),
            hs,
            lev + 1
          ));
          parent = newParent;
        } else {
          parent = n.children[childIndex];
        }
      } else {
        n.childrenBitmask[childIndex] = true;
      }
    }
  }

  draw(ctx: CanvasRenderingContext2D, nbs?: QuadTreeNode[]) {
    this.data[0][0].draw(ctx, 0, 0, ctx.canvas.width, 0, this.data, nbs);
  }
  
  getAsTravNode(x: number, y: number) {
    const rootTrav = new TravNode(null, 0, 0, 0, this.data);
    let lastTravNode = rootTrav;

    let parent = 0; // start with root
    let parentStartX = 0, parentStartY = 0;
    for (let lev = 0; lev < this.nLevels; lev++) {
      const n = this.data[lev][parent];
      const childCellSize = Math.pow(2, this.nLevels - lev - 1); // lev 0 -> 8, lev 1 -> 4, lev 2 -> 2, ...
      let childIndex = 0;
      if (x >= parentStartX + childCellSize) {
        childIndex += 1;
        parentStartX += childCellSize;
      }
      if (y >= parentStartY + childCellSize) {
        childIndex += 2;
        parentStartY += childCellSize;
      }
      if (n.childrenBitmask[childIndex]) {
        lastTravNode = new TravNode(lastTravNode, parent, childIndex, lev + 1, this.data);
        parent = n.children[childIndex];
      } else {
        break;
      }
    }
    return lastTravNode;
  }
}

////////////////////////
///// FLOOD FILL ///////
////////////////////////

class TravNode {
  constructor(
    public parent: TravNode,
    public parentIdx: number,
    public childIdx: number,
    public level: number,
    public readonly data: NodeData
  ) {

  }

  isLeaf() {
    if (this.level === this.data.length) return true;
    const node = this.getNode();
    return node === null || !node.hasChildren();
  }

  getNode() {
    if (this.level === 0) return this.data[0][0];
    const parent = this.getParentNode();
    if (parent.childrenBitmask[this.childIdx]) {
      return this.data[this.level][this.getParentNode().children[this.childIdx]];
    }
    return null;
  }

  getParentNode() {
    if (this.level === 0) throw new Error('Root has no parent');
    return this.data[this.level - 1][this.parentIdx];
  }

  getPos() {
    if (this.level === 0) {
      return [0, 0, ctx.canvas.width];
    } else {
      const pos = this.parent.getPos();
      const hs = pos[2] / 2;
      pos[0] += (this.childIdx % 2) * hs;
      pos[1] += Math.floor(this.childIdx / 2) * hs;
      pos[2] /= 2;
      return pos;
    }
  }

  draw(stroke: boolean, fill?: boolean) {
    const pos = this.getPos();
    // ctx.strokeStyle = "red";
    if (stroke) ctx.strokeRect(pos[0], pos[1], pos[2], pos[2]);
    if (fill) ctx.fillRect(pos[0], pos[1], pos[2], pos[2]);
  }
}

const getChildIndexForDirection = (dir: number, i: number) => {
  let childIdx = i * 2; // NX -> 0 or 2
  if (dir % 2 === 1) childIdx = i; // NY -> 0 or 1
  if (dir === 2) childIdx += 1; // PX -> 1 or 3
  if (dir === 3) childIdx += 2; // PY -> 2 or 3
  return childIdx % 4;
}

const getNeighbourGrtrEqSz = (tn: TravNode, dir: number, data: NodeData): TravNode | null => {
  if (tn.level === 0) return null;

  // For each of the 2 children on a side
  for (let i = 0; i < 2; ++i) {
    // If tn is a child at positive side, return a neighbour within the same parent at the opposite side, vice versa
    const oppositeSideChildIdx = getChildIndexForDirection(((dir + 2) % 4), i);
    // console.log('dir', dir, 'childIdx', tn.childIdx, 'i', i, 'opposite: ', oppositeSideChildIdx);
    if (tn.childIdx == oppositeSideChildIdx) {
      const curSideChildIdx = getChildIndexForDirection(dir, i);
      return new TravNode(tn.parent, tn.parentIdx, curSideChildIdx, tn.level, data);
    }
  }
  // Else, try to find the neighbour of the parent in this direction
  const pNb = getNeighbourGrtrEqSz(tn.parent, dir, data);
  // Return it if it's the nullNode, or it's a leaf
  if (pNb === null || pNb.isLeaf()) return pNb;

  const pNbIdx = data[pNb.level - 1][pNb.parentIdx].children[pNb.childIdx];
  // If the neighbour of the parent is not a leaf node, return its child that is the closest to the given node
  // For each of the 4 children on this side of the node...
  for (let i = 0; i < 2; ++i) {
    const curSideChildIdx = getChildIndexForDirection(dir, i);
    if (tn.childIdx == curSideChildIdx) {
      const oppositeSideChildIdx = getChildIndexForDirection(((dir + 2) % 4), i);
      return new TravNode(pNb, pNbIdx, oppositeSideChildIdx, pNb.level + 1, data);
    }
  }

  return null;
}

const getNeighboursSmSz = (tn: TravNode, tnNb: TravNode, dir: number, neighbours: TravNode[], data: NodeData) => {
  const candidates: TravNode[] = [];
  if (tnNb.parent != null) candidates.push(tnNb);

  while (candidates.length !== 0) {
    const can = candidates.pop();

    if (can.isLeaf()) {
      // if it's a leaf, add as a neighbour
      neighbours.push(can);
    } else {
      // else, add children in on the opposite side of the given direction as candidates
      const canIdx = can.getParentNode().children[can.childIdx];
      const newLev = can.level + 1;

      // For each of the 4 children on a side of the node, add candidates
      for (let i = 0; i < 2; ++i) {
        const oppositeSideChildIdx = getChildIndexForDirection(((dir + 2) % 4), i);
        candidates.push(new TravNode(can, canIdx, oppositeSideChildIdx, newLev, data));
      }
    }
  }
};

const getNeighbours = (tn: TravNode, neighbours: TravNode[], data: NodeData) => {
  // Get neighbours for all directions
  for (let dir = 0; dir < 4; ++dir) {
    // Get the neighbour directly next to it in this direction
    const nb = getNeighbourGrtrEqSz(tn, dir, data);
    // find the leaf nodes on the side of the input node inside the neighbour
    if (nb) {
      // neighbours.push(nb);
      getNeighboursSmSz(tn, nb, dir, neighbours, data);
    }
  }
};

const floodFill = (data: NodeData) => {
  const nLevels = data.length;
  ///// Finding node to start the floodfill /////
  const rootTrav = new TravNode(null, 0, 0, 0, data);
  let startTrav: TravNode | undefined;
  let prevParent = rootTrav;

  // Find the deepest node at the most negative corner of the scene
  const startIdx = 3; // Start flood fill at the top right of the scene
  let curNodeIdx = 0;
  for (let lev = 0; lev < nLevels; ++lev) {
    const node = data[lev][curNodeIdx];

    // Create new TravNode containing parent info
    const curNode = new TravNode(prevParent, curNodeIdx, startIdx, lev + 1, data);
    prevParent = curNode;
    curNodeIdx = data[lev][curNodeIdx].children[startIdx];
    // Replace the start node with its child
    startTrav = curNode;

    // Stop if no child is to be found
    if (!node.childrenBitmask[startIdx]) {
      break;
    }
  }

  // Check if the starting node intersects with geometry: Then the initial node is inside, so not a good starting point
  if (data[startTrav.level - 1][startTrav.parentIdx].childrenBitmask[startTrav.childIdx]) {
    throw new Error("Deepest corner node intersects with geometry, not dealing with this now... Aborting flood fill\n");
  }

  // Flood fill, using a queue of nodes that need to be visited
  const queue = [];
  queue.push(startTrav);
  data[startTrav.level - 1][startTrav.parentIdx][startTrav.childIdx] = true;

  let numOutsideLeafs = 0;

  while (queue.length !== 0) {
    const travNode = queue.pop();
    
    render();
    ctx.strokeStyle = "red";
    ctx.fillStyle = "purple";
    travNode.draw(false, true);


    data[travNode.level - 1][travNode.parentIdx].outsideMask[travNode.childIdx] = true;

    // For all neighbours
    const neighbours: TravNode[] = [];
    getNeighbours(travNode, neighbours, data);
    for (const nbTravNode of neighbours) {
      const nbP = nbTravNode.parent.getNode(); // neighbour parent
      // Node should be checked if 1. not intersects with geometry and 2. not already marked as outside
      if (!nbP.outsideMask[nbTravNode.childIdx]) {

        // TODO: This marks all neighbors as 'visited', so their neighbors are skipped when they are popped from the queue, right?

        // render(neighbours.map((tn) => tn.getNode()));
        // render();
        nbTravNode.draw(true, false);

        if (!nbP.childrenBitmask[nbTravNode.childIdx]) {
          // Add to queue to visit next
          queue.push(nbTravNode);
          numOutsideLeafs++;
        }
        else {
                  // Any neighbour of an outside node is also outside, as nodes on a surface also count as outside
          data[nbTravNode.level - 1][nbTravNode.parentIdx].outsideMask[nbTravNode.childIdx] = true;
        }
      }
    }
  }

  // Propagate it upwards - if all child is outside, set that outsideMask bit in the parent
  // for (let lev = nLevels - 2; lev >= 0; lev--) {
  //   for (let i = 0; i < data[lev].length; i++) {
  //     const n = data[lev][i];
  //     for (let c = 0; c < 4; c++) {
  //       if (n.childrenBitmask[c] && !data[lev + 1][n.children[c]].isInside()) {
  //         data[lev][i].outsideMask[c] = true;
  //       }
  //     }
  //   }
  // }
}


////////////////////////
////// MAIN ////////////
////////////////////////
const nLevels = 6;
const resolution = Math.pow(2, nLevels);

const canvas = document.querySelector<HTMLCanvasElement>('#canvas');
const ctx = canvas.getContext('2d');
ctx.lineWidth = 1;

const canvasSize = 1024;
canvas.width = canvasSize;
canvas.height = canvasSize;

const cellSize = canvasSize / resolution;

const tree = new QuadTree(nLevels);
let selectedNode: TravNode | undefined;

const render = (nbs?: QuadTreeNode[]) => {
  ctx.fillStyle = 'white';
  ctx.fillRect(0, 0, canvasSize, canvasSize);
  tree.draw(ctx, nbs);

  if (selectedNode) {
    selectedNode.getNode().drawDirect(ctx, 16, 'yellow')
    const nbs: TravNode[] = [];
    getNeighbours(selectedNode, nbs, tree.data);
    ctx.strokeStyle = 'blue';
    ctx.lineWidth = 6;
    nbs.forEach((nb) => nb.draw(true));
  }
}

const draw = (e: MouseEvent) => {
  if (!e.buttons) return;
  const x = Math.floor(e.offsetX / canvasSize * resolution);
  const y = Math.floor(e.offsetY / canvasSize * resolution);
  if (e.ctrlKey) {
    selectedNode = tree.getAsTravNode(x, y);
    console.log(selectedNode);
  } else {
    selectedNode = undefined;
    tree.setVoxel(x, y);
  }

  // ctx.fillStyle = 'black';
  // ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);

  render();
}

canvas.addEventListener('mousedown', draw);
canvas.addEventListener('mousemove', draw);

canvas.addEventListener('keydown', e => {
  console.log(e, e.key)
})

render();

// Tool box event handlers
const handleFloodFill = () => {
  floodFill(tree.data);
  render();
}

/**
 * TODO:
 * - Interactiveity
 * - - Ctrl click/right click to view node data
 * - - Option to render neighbors
 */
