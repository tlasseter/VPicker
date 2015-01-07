package com.Sts.Rescue;

/******************************************************************

  TestBed.java

  Testbed for RC squared Java to Rescue interface.  Note I am
  not making any attempt to generate a consistent model, just
  running some data thru the algorithms.

  Rod Hanks   April, 2000

*******************************************************************/
class TestBed
{
  public static void main(String[] args)
  {
//    System.out.println(RescueModel.BuildDate());
    RescueCoordinateSystem ultimateCs = new
            RescueCoordinateSystem("Austrailian TM 146 E",
                                    RescueCoordinateSystem.LDF, null,
                                     "northing", "m",
                                     "easting", "m",
                                     "elevation", "m");
    RescueVertex modelVertex = new RescueVertex("model vertex",
                                                ultimateCs, 96.5, 27.2, 1000);
    RescueCoordinateSystem penUltimateCs = new
                RescueCoordinateSystem("penultimate cs",
                  RescueCoordinateSystem.LDF, modelVertex,
                                     "northing", "m",
                                     "easting", "m",
                                     "elevation", "m");
		RescueModel model = new RescueModel("test model 1", penUltimateCs);

		model.ArchiveModel("tmodelx", false);
/**************************************************

  Note that in order to avoid holding the entire model in memory at once
  we have to first archive the model, then build the separately archivable
  parts, unload them, then archive again to make sure we got any related
  changes.

******************************************************/

    float xs[] = {0, 100, 100, 0, 25, 75, 75, 25};
    float ys[] = {0, 0, 100, 100, 25, 25, 75, 75};
    float zs[] = {0, 1, 2, 3, 4, 5, 6, 7};
    int howMany = xs.length;
    RescueTrimVertex[] vertices = new RescueTrimVertex[howMany];
    for (int loop = 0; loop < howMany; loop++)
    {
      vertices[loop] = new RescueTrimVertex(model, xs[loop], ys[loop], zs[loop]);
    }

    int lines[][] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
    howMany = lines.length;
    RescuePolyLine[] pLines = new RescuePolyLine[howMany];
    int ndx = 0;
    for (int loop = 0; loop < howMany; loop++)
    {
      pLines[loop] = new RescuePolyLine(model, vertices[lines[loop][0]],
                                               vertices[lines[loop][1]]);
      float x = (xs[lines[loop][0]] + xs[lines[loop][1]]) / 2;
      float y = (ys[lines[loop][0]] + ys[lines[loop][1]]) / 2;
      float z = (zs[lines[loop][0]] + zs[lines[loop][1]]) / 2;
      RescuePolyLineNode middleNode = new RescuePolyLineNode(x, y, z);
      pLines[loop].AddPolyLineNode(middleNode);
    }

    RescueTrimEdge[] edges = new RescueTrimEdge[pLines.length];
    for (int loop = 0; loop < pLines.length; loop++)
    {
      edges[loop] = new RescueTrimEdge(pLines[loop], RescueTrimEdge.R_RIGHT_TO_LEFT);
    }

    RescueTrimLoop boundaryLoop = new RescueTrimLoop();
    for (int loop = 0; loop < 4; loop++)
    {
      boundaryLoop.AddLoopEdge(edges[loop]);
    }
    RescueTrimLoop interiorLoop = new RescueTrimLoop();
    for (int loop = 0; loop < 4; loop++)
    {
      interiorLoop.AddLoopEdge(edges[loop + 4]);
    }

    RescueSection section = new RescueSection(RescueCoordinateSystem.LUB,
											                        "test section",
											                        model,
											                        RescueSurface.AUXILLIARY,
											                        0, 4, 0, 4, (float) -999.25);
    RescueEdgeSet edgeSet = section.Edges();
    edgeSet.AddBoundaryLoop(boundaryLoop);
    edgeSet.AddInteriorLoop(interiorLoop);

    for (int loop = 0; loop < vertices.length; loop++)
    {
      vertices[loop].SetUVValue(section, (float) (0.5 * (float) loop),
                                         (float) (2.5 * (float) loop));
    }

    RescueTripletArray geometry = section.Geometry();
    float[] sx = new float[16];
    float[] sy = new float[16];
    float[] sz = new float[16];
    for (int loop = 0; loop < 16; loop++)
    {
      sx[loop] = (float) (0.5 * (float) loop);
      sy[loop] = (float) (100 * (float) loop);
      sz[loop] = (float) (10000 * (float) loop);
    }
    geometry.AssignXValue(sx);
    geometry.AssignYValue(sy);
    geometry.AssignZValue(sz);
    geometry.Unload();

    RescueHorizon topH = new RescueHorizon("Top Horizon", model);
    RescueHorizon bottomH = new RescueHorizon("Bottom Horizon", model);
    RescueUnit unit = new RescueUnit("Unit 1", model);
    RescueBlock block = new RescueBlock("Block 1", model);
    topH.SetUnitBelowMe(unit);
    bottomH.SetUnitAboveMe(unit);
    unit.SetHorizonAboveMe(topH);
    unit.SetHorizonBelowMe(bottomH);
    RescueBlockUnitHorizonSurface topSurf = new RescueBlockUnitHorizonSurface(
                                         RescueCoordinateSystem.LUB,
                                         topH, 19000, 10, 0, 4,
                                         22000, (float) 12.5, 0, 4,
                                         (float) -999.25, 0, null,
                                         RescueSurface.HORIZON);
    RescueBlockUnitHorizonSurface bottomSurf = new RescueBlockUnitHorizonSurface(
                                         RescueCoordinateSystem.LUB,
                                         bottomH, 19000, 10, 0, 4,
                                         22000, (float) 12.5, 0, 4,
                                         (float) -999.25, 0, null,
                                         RescueSurface.HORIZON);
/*
    RescueBlockUnit bu = new RescueBlockUnit(RescueCoordinateSystem.LUB,
                                              block, unit, 19000, 10, 0, 4,
                                              22000, (float) 12.5, 0, 4, 0, 4,
                                              (float) -999.25, 0, null);
*/
/*
    RescueBlockUnit bu = new RescueBlockUnit(RescueCoordinateSystem.LUB,
                                              block, unit, 19000, 10, 0, 4,
                                              22000, (float) 12.5, 0, 4, 0, 4,
                                              (float) -999.25,
                                              topSurf, 0, bottomSurf, 0,
                                              0, null);
*/
    RescueBlockUnit bu = new RescueBlockUnit(RescueCoordinateSystem.LUB,
                                              block, unit, 19000, 10, 0, 4,
                                              22000, (float) 12.5, 0, 4, 0, 4,
                                              (float) -999.25,
                                              topSurf, 0, 10,
                                              RescueGeometry.R_ONLAP,
                                              0, null);
    topSurf.SetBlockUnitAboveMe(bu);
    bottomSurf.SetBlockUnitBelowMe(bu);

    RescueTripletArray topGeom = topSurf.Geometry();
    for (int loop = 0; loop < 16; loop++)
    {
      sz[loop] = (float) (5000 + (float) loop);
    }
    topGeom.AssignZValue(sz);
    topGeom.Unload();

    RescueTripletArray bottomGeom = bottomSurf.Geometry();
    for (int loop = 0; loop < 16; loop++)
    {
      sz[loop] = (float) (7000 + (float) loop);
    }
    bottomGeom.AssignZValue(sz);
    bottomGeom.Unload();

    RescueBlockUnitPropertyGroup group = new RescueBlockUnitPropertyGroup("myGroup", bu);

    float[][][] propArray = new float[3][3][3];
    for (int layer = 0; layer < 3; layer++)
    {
      for (int row = 0; row < 3; row++)
      {
        for (int col = 0; col < 3; col++)
        {
          propArray[layer][row][col] = (layer * 10000) + (row * 1000) + (col * 100);
        }
      }
    }
    RescueProperty property = new RescueProperty(group, "Kx", "porosity", "decimal",
                                                 (float) -999.25, propArray);
    property.Data().Unload();

    RescueMacroVolume mvol = bu.MacroVolume();
    RescueSection sideSection = new RescueSection(RescueCoordinateSystem.LUB,
											                        "side section",
											                        model,
											                        RescueSurface.AUXILLIARY,
											                        0, 4, 0, 4, (float) -999.25);
    RescueBlockUnitSide unitSide = new RescueBlockUnitSide(sideSection);
    mvol.AddBlockUnitSide(unitSide);
    RescueSection insideSection = new RescueSection(RescueCoordinateSystem.LUB,
											                        "inside section",
											                        model,
											                        RescueSurface.AUXILLIARY,
											                        0, 4, 0, 4, (float) -999.25);
    mvol.AddInteriorSection(insideSection);
    RescueEdgeSet sideEdges = unitSide.Edges();
    GenerateEdgeSet(model, sideEdges, 100);

    RescueEdgeSet topEdges = new RescueEdgeSet();
    GenerateEdgeSet(model, topEdges, 200);
    mvol.SetTopEdge(topEdges);

    RescueEdgeSet bottomEdges = new RescueEdgeSet();
    GenerateEdgeSet(model, bottomEdges, 300);
    mvol.SetBottomEdge(topEdges);

    for (int layer = 0; layer < 4; layer++)
    {
      RescueEdgeSet layerEdge = new RescueEdgeSet();
      GenerateEdgeSet(model, layerEdge, (layer + 1) * 1000);
      mvol.AddKLayerEdge(layerEdge);
    }

    model.UnloadWireframe();
    model.ArchiveModel();
    model.dispose();
  }

  private static void GenerateEdgeSet(RescueModel model,
                                      RescueEdgeSet edgeSet,
                                      float data)
  {
    float xs[] = {0, 100, 100, 0, 25, 75, 75, 25};
    float ys[] = {0, 0, 100, 100, 25, 25, 75, 75};
    float zs[] = {0, 1, 2, 3, 4, 5, 6, 7};
    int howMany = xs.length;
    RescueTrimVertex[] vertices = new RescueTrimVertex[howMany];
    for (int loop = 0; loop < howMany; loop++)
    {
      vertices[loop] = new RescueTrimVertex(model, xs[loop], ys[loop], zs[loop] + data);
    }

    int lines[][] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
    howMany = lines.length;
    RescuePolyLine[] pLines = new RescuePolyLine[howMany];
    int ndx = 0;
    for (int loop = 0; loop < howMany; loop++)
    {
      pLines[loop] = new RescuePolyLine(model, vertices[lines[loop][0]],
                                               vertices[lines[loop][1]]);
      float x = (xs[lines[loop][0]] + xs[lines[loop][1]]) / 2;
      float y = (ys[lines[loop][0]] + ys[lines[loop][1]]) / 2;
      float z = (zs[lines[loop][0]] + zs[lines[loop][1]]) / 2;
      RescuePolyLineNode middleNode = new RescuePolyLineNode(x, y, z);
      pLines[loop].AddPolyLineNode(middleNode);
    }

    RescueTrimEdge[] edges = new RescueTrimEdge[pLines.length];
    for (int loop = 0; loop < pLines.length; loop++)
    {
      edges[loop] = new RescueTrimEdge(pLines[loop], RescueTrimEdge.R_LEFT_TO_RIGHT);
    }

    RescueTrimLoop boundaryLoop = new RescueTrimLoop();
    for (int loop = 0; loop < 4; loop++)
    {
      boundaryLoop.AddLoopEdge(edges[loop]);
    }
    RescueTrimLoop interiorLoop = new RescueTrimLoop();
    for (int loop = 0; loop < 4; loop++)
    {
      interiorLoop.AddLoopEdge(edges[loop + 4]);
    }
    edgeSet.AddBoundaryLoop(boundaryLoop);
    edgeSet.AddInteriorLoop(interiorLoop);
  }
}
