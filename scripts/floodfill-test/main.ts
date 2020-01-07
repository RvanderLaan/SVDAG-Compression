
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

  constructor() {
    this.childrenBitmask = []
    this.children = [];
    this.outsideMask = [];
    for (let i = 0; i < NCHILDREN; i++) {
      this.childrenBitmask[i] = false;
      this.children[i] = nullNode;
      this.outsideMask[i] = false;
    }
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

  async draw(ctx: CanvasRenderingContext2D, x: number, y: number, size: number, lev: number, data: NodeData) {
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

    if (this.isInside()) {
      ctx.fillStyle = 'rgba(0,0,200,0.5)';
      ctx.fillRect(x, y, size, size);
    }
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
    this.data[0].push(new QuadTreeNode()); // root node
  }

  setVoxel(x: number, y: number) {
    const r = Math.pow(2, this.nLevels);
    let parent = 0; // start with root
    let parentStartX = 0, parentStartY = 0;
    for (let lev = 0; lev < this.nLevels; lev++) {
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
        if (!this.data[lev][parent].childrenBitmask[childIndex]) {
          const newParent = this.data[lev + 1].length;
          this.data[lev][parent].setChild(childIndex, newParent);
          this.data[lev + 1].push(new QuadTreeNode());
          parent = newParent;
        } else {
          parent = this.data[lev][parent].children[childIndex];
        }
      } else {
        this.data[lev][parent].childrenBitmask[childIndex] = true;
      }
    }
  }

  draw(ctx: CanvasRenderingContext2D) {
    this.data[0][0].draw(ctx, 0, 0, ctx.canvas.width, 0, this.data);
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
    return node !== null && !node.hasChildren();
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
}

const getChildIndexForDirection = (dir: number, i: number) => {
  let childIdx = i; // NX
  if (dir % 2 === 1) childIdx = i + 1; // NY
  if (dir === 2) childIdx += 2; // PX
  if (dir === 3) childIdx += 1; // PY
  return childIdx;
}

const getNeighbourGrtrEqSz = (tn: TravNode, dir: number, data: NodeData): TravNode | null => {
  if (tn.level === 0) return null;

  // For each of the 2 children on a side
  for (let i = 0; i < 2; ++i) {
    // If tn is a child at positive side, return a neighbour within the same parent at the opposite side, vice versa
    const oppositeSideChildIdx = getChildIndexForDirection(((dir + 2) % 4), i);
    if (tn.childIdx == oppositeSideChildIdx) {
      const curSideChildIdx = getChildIndexForDirection(dir, i);
      return new TravNode(tn.parent, tn.parentIdx, curSideChildIdx, tn.level, data);
    }
  }
  // Else, try to find the neighbour of the parent in this direction
  const pNb = getNeighbourGrtrEqSz(tn.parent, dir, data);
  // Return it if it's the nullNode, or it's a leaf
  if (pNb.parent === null || pNb.isLeaf()) return pNb;

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

const candidates: TravNode[] = [];
const getNeighboursSmSz = (tn: TravNode, tnNb: TravNode, dir: number, neighbours: TravNode[], data: NodeData) => {
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
    getNeighboursSmSz(tn, nb, dir, neighbours, data);
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

    console.log(queue);

    data[travNode.level - 1][travNode.parentIdx].outsideMask[travNode.childIdx] = true;

    // For all neighbours
    const neighbours: TravNode[] = [];
    getNeighbours(travNode, neighbours, data);
    for (const nbTravNode of neighbours) {
      const nbP = nbTravNode.parent.getNode(); // neighbour parent
      // Node should be checked if 1. not intersects with geometry and 2. not already marked as outside
      if (!nbP.outsideMask[nbTravNode.childIdx]) {

        // This marks all neighbors as 'visited', so their neighbors are skipped when they are popped from the queue, right?

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

const render = () => {
  ctx.fillStyle = 'white';
  ctx.fillRect(0, 0, canvasSize, canvasSize);
  tree.draw(ctx);
}



const draw = (e: MouseEvent) => {
  if (!e.buttons) return;
  const x = Math.floor(e.offsetX / canvasSize * resolution),
    y = Math.floor(e.offsetY / canvasSize * resolution);

  // ctx.fillStyle = 'black';
  // ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);

  tree.setVoxel(x, y);
  render();
}

canvas.addEventListener('mousedown', draw);
canvas.addEventListener('mousemove', draw);

canvas.addEventListener('keydown', e => {
  console.log(e, e.key)
})

render();

// Tool box event handlers
const handleFloodFill = () => floodFill(tree.data);
