
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;

public class StsBlockSide extends StsObject
{
    StsSection section;
    int side;
    StsLine prevLine;
    StsLine nextLine;
    StsBlock block; // assigned when blockSide is added to block
    
    transient StsList surfaceEdges;

    static final int LEFT = StsParameters.LEFT;
    static final int RIGHT = StsParameters.RIGHT;
    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    public StsBlockSide()
         {
         }

    public StsBlockSide(StsLine prevLine, StsLine nextLine, StsSection section, int side)
         {
             this.prevLine = prevLine;
             this.nextLine = nextLine;
             this.section = section;
             this.side = side;
         }

    public boolean initialize(StsModel model) { return true; }

    public StsLine getPrevLine() { return prevLine; }
    public StsLine getNextLine() { return nextLine; }
    public void setBlock(StsBlock block)
         {
             fieldChanged("block", block);
         }

    public StsBlock getBlock() { return block; }
    public StsSection getSection() { return section; }
    public int getSide() { return side; }

    public StsBlockSide getNextBlockSide()
         {
             StsSection onSection = nextLine.getOnSection();
             if(onSection != null)
             {
                 if(onSection == section)  // right turn off of section onto connected section
                 {
                     StsSection connectedSection = nextLine.getOnlyConnectedSection();
                     if(connectedSection == null) return null;
                     return connectedSection.getNextBlockSide(nextLine);
                 }
                 else // right turn off of section onto onSection
                     return onSection.getNextBlockSide(nextLine);
             }
             StsSection nextSection = nextLine.getOtherSection(section);
             if(nextSection == null)
             {
                 StsException.systemError("StsBlockSide.getNextSide() failed."+
                     " Couldn't find next blockSide from blockSide: " + getLabel());
                 return null;
             }
             return nextSection.getNextBlockSide(nextLine);
         }


    public StsBlockSide getPrevBlockSide()
         {
             StsSection onSection = prevLine.getOnSection();
             if(onSection != null)
             {
                 if(onSection == section)  // left turn off of section onto connected section
                 {
                     StsSection connectedSection = prevLine.getOnlyConnectedSection();
                     if(connectedSection == null) return null;
                     return connectedSection.getPrevBlockSide(prevLine);
                 }
                 else // left turn off of section onto onSection
                     return onSection.getPrevBlockSide(prevLine);
             }
             StsSection prevSection = prevLine.getOtherSection(section);
     //        StsSection prevSection = prevLine.getPrevSection();
             if(prevSection == null)
             {
                 StsException.systemError("StsBlockSide.getPrevSide() failed."+
                     " Couldn't find prev blockSide from blockSide: " + getLabel());
                 return null;
             }
             return prevSection.getPrevBlockSide(prevLine);
         }

    public StsSurfaceEdge constructInitialEdge(StsModelSurface modelSurface, StsBlockGrid blockGrid)
         {
             StsSurfaceVertex vertex = prevLine.getConstructSurfaceEdgeVertex(modelSurface, block);
             if(vertex == null) return null;

             StsSurfaceVertex nextVertex = nextLine.getConstructSurfaceEdgeVertex(modelSurface, block);
             if(nextVertex == null) return null;

             StsSurfaceEdge edge = getSurfaceEdge(vertex, nextVertex);
             if(edge == null)
             {
                 edge = new StsSurfaceEdge(vertex, nextVertex, modelSurface, section, side, blockGrid);
                 addEdge(edge);
             }
             else
                 edge.setBlockGrid(blockGrid);
             edge.constructSurfaceEdgePoints(modelSurface, false);
             return edge;
         }

    private StsSurfaceEdge getSurfaceEdge(StsSurfaceVertex prevVertex, StsSurfaceVertex nextVertex)
         {
             if(surfaceEdges == null) return null;
             int nSurfaceEdges = surfaceEdges.getSize();
             for(int n = 0; n < nSurfaceEdges; n++)
             {
                 StsSurfaceEdge edge = (StsSurfaceEdge)surfaceEdges.getElement(n);
                 if(edge.prevVertex == prevVertex && edge.nextVertex == nextVertex)
                     return edge;
             }
             return null;
         }

    public void addEdge(StsSurfaceEdge surfaceEdge)
         {
             if(surfaceEdges == null) surfaceEdges = new StsList(4, 2);
             surfaceEdges.add(surfaceEdge);
         }

    public StsSurfaceEdge getSurfaceEdge(StsModelSurface surface)
         {
             if(surfaceEdges == null) return null;
             for(int n = 0; n < surfaceEdges.getSize(); n++)
             {
                 StsSurfaceEdge edge = (StsSurfaceEdge)surfaceEdges.getElement(n);
                 if(edge.getSurface() == surface) return edge;
             }
             return null;
         }

    public String getLabel()
         {
             return new String("BlockSide on section: " + section.getLabel() + " side: " +
                 StsParameters.sideLabel(side) + " from: " + prevLine.getLabel() + " to: " +
                 nextLine.getLabel());
         }

    public String toString()
         {
             return new String("BlockSide on section: " + section.getLabel() + " side: " +
                 StsParameters.sideLabel(side) + " from: " + prevLine.getLabel() + " to: " +
                 nextLine.getLabel());
         }

    public boolean delete()
    {
        super.delete();
        surfaceEdges = null;
        prevLine.deleteSurfaceVertices();
        nextLine.deleteSurfaceVertices();
        return true;
    }
}
