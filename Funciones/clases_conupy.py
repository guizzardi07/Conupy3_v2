from collections.abc import MutableSequence
import numpy as np
from Funciones.script_funciones import dist

class NodesList(MutableSequence):

    def __init__(self):
        self.list = list()
        self.geo_hash = {}

    def __len__(self):
        return len(self.list)

    def __getitem__(self, i):
        return self.list[i]

    def __delitem__(self, i):
        del self.list[i]

        # Reconstruct the geo_hash
        self.geo_hash = {}
        for i, v in enumerate(self.list):
            ix, iy = [int(z/100.0) for z in v.p]
            self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
            if not i in self.geo_hash[(ix,iy)]:
                self.geo_hash[(ix,iy)].append(i)

    def __setitem__(self, i, v):
        # Update geo cache
        ix, iy = [int(z/100.0) for z in v.p]
        self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
        if not i in self.geo_hash[(ix,iy)]:
            self.geo_hash[(ix,iy)].append(i)

        self.list[i] = v

    def insert(self, i, v):
        # TODO: don't allow insertions in between

        # Update geo cache
        ix, iy = [int(z/100.0) for z in v.p]
        self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
        if not i in self.geo_hash[(ix,iy)]:
            self.geo_hash[(ix,iy)].append(i)

        self.list.insert(i, v)

    def __str__(self):
        return str(self.list)

    def getINodesInside(self, p0, p1):
        ix0, iy0 = [int(z/100.0) for z in p0]
        ix1, iy1 = [int(z/100.0) for z in p1]
        inodes = []
        for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
            for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                if (ix, iy) in self.geo_hash:
                    inodes.extend(self.geo_hash[(ix, iy)])
        return inodes

    def getINodesNear(self, p0, distance):
        return self.getINodesInside(np.add(p0, -distance), np.add(p0, distance))

    def getINode(self, p0):
        inodes = self.getINodesInside(p0, p0)
        mindist, minj = 0.5, -1
        for j in inodes:
            nodo2 = self.list[j]
            d = dist(p0, nodo2.p)
            if d < mindist:
                mindist = d
                minj = j
        if minj == -1:
            return None
        return minj

class LinksDict(dict):
        def __init__(self, nodes):
            self.nodes = nodes
            self.geo_hash = {}
            super(dict, self).__init__()

        def __setitem__(self, key, value):
            # optional processing here

            # Update geo cache
            n0, n1 = key
            ix0, iy0 = [int(z/100.0) for z in self.nodes[n0].p]
            ix1, iy1 = [int(z/100.0) for z in self.nodes[n1].p]
            for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
                for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                    self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
                    self.geo_hash[(ix,iy)].append((n0, n1, value))

            super(LinksDict, self).__setitem__(key, value)

        def __delitem__(self, key):
            # Update geo cache
            value = self[key]
            n0, n1 = key
            ix0, iy0 = [int(z/100.0) for z in self.nodes[n0].p]
            ix1, iy1 = [int(z/100.0) for z in self.nodes[n1].p]
            for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
                for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                    self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
                    self.geo_hash[(ix,iy)].remove((n0, n1, value))

            return super(LinksDict, self).__delitem__(key)

        def getLinksInside(self, p0, p1):
            ix0, iy0 = [int(z/100.0) for z in p0]
            ix1, iy1 = [int(z/100.0) for z in p1]
            links = dict()
            for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
                for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                    if (ix, iy) in self.geo_hash:
                        for n0, n1, link in self.geo_hash[(ix, iy)]:
                            links[(n0, n1)] = link
            return links

        def getLinksNear(self, p0, distance):
            return self.getLinksInside(np.add(p0, -distance), np.add(p0, distance))

        def update(self, *args, **kwargs):
            if args:
                if len(args) > 1:
                    raise TypeError("update expected at most 1 arguments, "
                                    "got %d" % len(args))
                other = dict(args[0])
                for key in other:
                    self[key] = other[key]
            for key in kwargs:
                self[key] = kwargs[key]

        def setdefault(self, key, value=None):
            if key not in self:
                self[key] = value
            return self[key]
