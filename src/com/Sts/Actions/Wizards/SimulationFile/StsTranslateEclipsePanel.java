

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.Simulation.*;
import com.Sts.MVC.*;

import java.util.*;


public class StsTranslateEclipsePanel extends StsJPanel
{
    private int eclRow, eclCol, eclLayer;
    private int row, col, layer, nBlock;
    private float x, y, z;
    StsBlock[] blocks;
    int nBlocks;
    int nEclipseRows;

    public StsTranslateEclipsePanel(StsModel model)
    {
        try
        {
            initializeEclipseParameters(model);
            buildPanel();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void initializeEclipseParameters(StsModel model)
    {
        StsBuiltModel builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
        blocks = builtModel.getBlocks();
        nBlocks = blocks.length;
        StsBlock lastBlock = blocks[nBlocks-1];
        nEclipseRows = lastBlock.eclipseRowMin + lastBlock.nCellRows - 1;
    }

    private void buildPanel() throws Exception
    {
        gbc.fill = HORIZONTAL;

        StsGroupBox eclipseIndexBox = new StsGroupBox("Eclipse Indexes");

        StsIntFieldBean eclRowBean = new StsIntFieldBean(this, "eclRow", 1, 1000, "Eclipse Row J");
        eclipseIndexBox.addBeanToRow(eclRowBean);
        StsIntFieldBean eclColBean = new StsIntFieldBean(this, "eclCol", 1, 1000, "Eclipse Col I");
        eclipseIndexBox.addBeanToRow(eclColBean);
        StsIntFieldBean eclLayerBean = new StsIntFieldBean(this, "eclLayer", 1, 1000, "Eclipse Layer K");
        eclipseIndexBox.addBeanEndRow(eclLayerBean);

        add(eclipseIndexBox);

        StsGroupBox s2sIndexBox = new StsGroupBox("S2S Indexes");
        StsIntFieldBean rowBean = new StsIntFieldBean(this, "row", -1, 1000, "S2S Row");
        s2sIndexBox.addBeanToRow(rowBean);
        StsIntFieldBean colBean = new StsIntFieldBean(this, "col", -1, 1000, "S2S Col");
        s2sIndexBox.addBeanToRow(colBean);
        StsIntFieldBean layerBean = new StsIntFieldBean(this, "layer", -1, 1000, "S2S Layer");
        s2sIndexBox.addBeanToRow(layerBean);
        StsIntFieldBean blockBean = new StsIntFieldBean(this, "nBlock", -1, 1000, "S2S Block");
        s2sIndexBox.addBeanEndRow(blockBean);

        add(s2sIndexBox);

        StsButton cellPropertiesButton = new StsButton("Cell Summary", "Provides info on cell geometry.", this, "cellProperties");
        gbc.fill = NONE;
        add(cellPropertiesButton);
    }

    public void cellProperties()
    {
        StsBlock block = getBlock();
        if(block == null) return;
        StsBlock.BlockCellColumn.GridCell gridCell = block.getParentCellOrGridCell(row, col, layer);
        if(gridCell == null) return;

        StsTextAreaDialog textAreaDialog = new StsTextAreaDialog(null, "Grid Cell Properties", false);
        textAreaDialog.appendLine(gridCell.toString());
        ArrayList<StsNNC> nncList = gridCell.nncList;
        if(nncList != null)
        {
            textAreaDialog.appendLine("");
            for(StsNNC nnc : nncList)
                textAreaDialog.appendLine(nnc.toString());
        }
        textAreaDialog.setVisible(true);
    }

    private StsBlock getBlock()
    {
        if(nBlock < 0 || nBlock >= nBlocks) return null;
        return blocks[nBlock];
    }

    public void setEclRow(int eclRow)
    {
        this.eclRow = eclRow;
        computeS2sIndexes();
    }

    public void setEclCol(int eclCol)
    {
        this.eclCol = eclCol;
        computeS2sIndexes();
    }

    public void setEclLayer(int eclLayer)
    {
        this.eclLayer = eclLayer;
        computeS2sIndexes();
    }

    private void computeS2sIndexes()
    {
        for(int n = 0; n < nBlocks; n++)
        {
            StsBlock block = blocks[n];
            int eclRowMin = block.eclipseRowMin;
            int eclRowMax = eclRowMin + block.nCellRows - 1;
            if(eclRow >= eclRowMin && eclRow <= eclRowMax)
            {
                row = block.rowMax - (eclRow - eclRowMin) - 1;
                col = eclCol + block.colMin - 1;
                layer = eclLayer - 1;
                nBlock = n;
                updateBeans();
                return;
            }
        }
        row = -1;
        col = -1;
        layer = eclLayer - 1;
        nBlock = -1;
        updateBeans();
    }
    public void setRow(int row)
    {
        this.row = row;
        computeS2sIndexes();
    }

    public void setCol(int col)
    {
        this.col = col;
        computeS2sIndexes();
    }

    public void setLayer(int layer)
    {
        this.layer = layer;
        computeEclIndexes();
    }

    public void setNBlock(int nBlock)
    {
        this.nBlock = nBlock;
        computeEclIndexes();
    }

    private void computeEclIndexes()
    {
        if(nBlock < 0 || nBlock >= nBlocks)
        {
            eclRow = -1;
            eclCol = -1;
            eclLayer = eclLayer + 1;
        }
        else
        {
            StsBlock block = blocks[nBlock];
            int rowMin = block.rowMin;
            int rowMax = block.rowMax - 1;
            if(row >= rowMin && row <= rowMax)
            {
                eclRow = rowMax - row + block.eclipseRowMin;
                eclCol = col + block.colMin + 1;
                eclLayer = eclLayer + 1;
            }
            else
            {
                eclRow = -1;
                eclCol = -1;
                eclLayer = eclLayer + 1;
            }
        }
        updateBeans();
    }

    public static void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsTranslateEclipsePanel translatePanel = new StsTranslateEclipsePanel(model);
        StsToolkit.createDialog("Translate Eclipse Test", translatePanel, true);
    }

    public int getEclRow()
    {
        return eclRow;
    }

    public int getEclCol()
    {
        return eclCol;
    }

    public int getEclLayer()
    {
        return eclLayer;
    }

    public int getRow()
    {
        return row;
    }

    public int getCol()
    {
        return col;
    }

    public int getLayer()
    {
        return layer;
    }

    public int getNBlock()
    {
        return nBlock;
    }
}