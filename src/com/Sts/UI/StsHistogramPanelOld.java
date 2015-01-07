//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;

public class StsHistogramPanelOld extends JPanel implements ActionListener
{
    public static final int VERTICAL = 0;
    public static final int HORIZONTAL = 1;
    int orientation = VERTICAL;

    ActionListener[] actionListeners = null;

    float verticalScale = -1.0f;
    JPanel histogramPanel = new JPanel();
    DecimalFormat colorFormat = new DecimalFormat("#.0##");
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    Insets insets = null;
    Font defaultFont = new Font("Dialog", 0, 11);
    Font floatFont = new Font("Dialog", 0, 9);
    FontMetrics fm = getFontMetrics(floatFont);
    int dataX = 0, dataY = 0, dataXOld = -999, dataYOld = -999, dataWidth = 40, dataHeight = 20;
    Image offScreenImage = null;
    Graphics offScreenGraphics;
    Float dataValue = null;
    Float dataPercent = null;
    float[] opacityValues = null;
    int[] opacityX = new int[257];
    int[] opacityY = new int[257];
    StsColorscale colorscale = null;
    StsColorscalePanel colorscalepanel = null;

    private JMenuItem kc = new JMenuItem("Set to Fully Opaque");

    int mouseBtn = 0;
    int previousIdx = -1;

    int maxWidth = 50;
    int nPoints = 255;
    int heightPerInterval = 1;
    int barHeight = 255;
    int barWidth = 50;
    Rectangle barRect;
    int xOrigin, yOrigin;
    float dataMax = -1000.0f;
    float dataMin = 1000.0f;
    float minVal = 1;
    float maxVal = 255;
    int minIdx = 0;
    int maxIdx = 254;
    float topClip = 0.0f;
    float btmClip = 0.0f;

    int numSegments = 4;
    boolean segmentsOn = false;

    float[] dataHist = null;

    public StsHistogramPanelOld()
    {
        this(VERTICAL, null);
    }

    public StsHistogramPanelOld(int orient)
    {
        this(orient, null);
    }

    public StsHistogramPanelOld(float[] histData)
    {
        this(VERTICAL, histData);
    }

    public StsHistogramPanelOld(int orient, float[] histData)
    {
        orientation = orient;
        try
        {
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
        if(histData == null)
            dataHist = new float[255];
        else
            dataHist = histData;
        init();
    }

    public void setData(float[] histData)
    {
        dataHist = histData;
        if (dataHist == null) return;
        init();
    }

    public void setVolRange(float min, float max)
    {
        minVal = min;
        maxVal = max;
    }

    public void setColorscalePanel(StsColorscalePanel colorscalepanel)
    {
        this.colorscale = colorscalepanel.getColorscale();
        this.colorscalepanel = colorscalepanel;
        resetOpacityValues();
    }

    private void init()
    {
        addMouseListener();
        addMouseMotionListener();
        if (dataHist.length != 255)
        {
            System.err.println("Data Histogram must have 255 elements");
        }
        for (int i = 0; i < dataHist.length; i++)
        {
            if (dataHist[i] > dataMax)
            {
                dataMax = dataHist[i];
            }
            if (dataHist[i] < dataMin)
            {
                dataMin = dataHist[i];
            }
        }
        resetOpacityValues();
        kc.addActionListener(this);
    }

    public void resetOpacityValues()
    {
        barRect = histogramPanel.getBounds();

        opacityX[0] = barRect.x + 2;
        opacityY[0] = nPoints + barRect.y;
        barRect = histogramPanel.getBounds();
        int bwidth = barRect.width;
        if(bwidth == 0)
            bwidth = barWidth;

		if(colorscale == null)
		{
			if (opacityValues == null) opacityValues = new float[nPoints + 1];
			for (int n = 0; n < nPoints; n++)
				opacityValues[n] = 1.0f;
		}
		else
        {
            int[] colorscaleOpacityValues = colorscale.getOpacityValues();
            if (opacityValues == null) opacityValues = new float[nPoints + 1];
            for(int n = 0; n < nPoints; n++)
                opacityValues[n] = colorscaleOpacityValues[n]/255.0f;
        }

        for (int i = 1; i < nPoints + 1; i++)
		{
			opacityX[i] = (int)(bwidth * opacityValues[i-1]) + 2;
            opacityY[i] = barRect.y + nPoints - i + 1;
        }
        opacityX[nPoints + 1] = barRect.x + 2;
        opacityY[nPoints + 1] = barRect.y;
        repaint();
    }

    public void clearOpacityValues()
    {
        barRect = histogramPanel.getBounds();

        opacityX[0] = barRect.x + 1;
        opacityY[0] = nPoints + barRect.y;
        barRect = histogramPanel.getBounds();
        int bwidth = barRect.width;
        if(bwidth == 0)
            bwidth = barWidth;

        for (int i = colorscalepanel.minSliderIndex + 1; i < colorscalepanel.maxSliderIndex + 1; i++)
        {
            if(colorscale == null)
                opacityValues[i - 1] = 1.0f;
            else
            {
                Color[] colors = colorscale.getNewColors();
                opacityValues[i - 1] = 1.0f;
            }
            opacityX[i] = (int)(bwidth * opacityValues[i-1]) + 1;
            opacityY[i] = barRect.y + nPoints - i + 1;
        }
        opacityX[nPoints + 1] = barRect.x + 1;
        opacityY[nPoints + 1] = barRect.y;
        repaint();

        if (colorscale != null)
            colorscalepanel.setOpacityValues(opacityValues);
    }

    private void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
        histogramPanel.setFont(null);
        Dimension histogramSize = new Dimension(barWidth, barHeight);
        if (orientation == HORIZONTAL)
        {
            histogramSize = new Dimension(barHeight, barWidth);

        }
        histogramPanel.setMinimumSize(histogramSize);
        histogramPanel.setPreferredSize(histogramSize);
        Dimension size = new Dimension(histogramSize.width + 4,
                                       histogramSize.height + 4);

        setMinimumSize(size);
        setPreferredSize(size);

        this.add(histogramPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE,
            new Insets(2, 2, 2, 2), 0, 0));
    }

    public void setVerticalScale(float scale)
    {
        verticalScale = scale;
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();
        // Process Radio Button Input
        if(source == kc)
            if(colorscale != null) clearOpacityValues();
    }

    public void paint(Graphics g)
    {
        if (dataHist == null)
        {
            return;
        }
        if (g == null)
        {
            return;
        }
        super.paint(g);

        barRect = histogramPanel.getBounds();
        xOrigin = barRect.x;
        yOrigin = barRect.y - 1;

        Rectangle back = this.getVisibleRect();
        g.setColor(getBackground());
        g.fillRect(back.x, back.y, back.width, back.height);
        paintChildren(g);

        int bwidth = barRect.width;
        if (offScreenImage == null)
        {
            dataWidth = fm.stringWidth("99999.99");
            dataWidth += dataWidth * 0.05;
            if (dataWidth > (bwidth - 2))
            {
                dataWidth = bwidth - 2;
            }
//            dataHeight = fm.getHeight() * 2;
            dataHeight = fm.getHeight();
            dataHeight += dataHeight * 0.05;
            offScreenImage = createImage(dataWidth, dataHeight);
        }
        offScreenGraphics = offScreenImage.getGraphics();

        // Draw the histogram
        int x = xOrigin;
        int y = yOrigin;
        if (orientation == HORIZONTAL)
        {
            x = xOrigin;
            y = yOrigin;
            bwidth = barRect.height;
        }

        // Draw Color Scale - Bottom up
        float scale = bwidth / dataMax;
        topClip = 0.0f;
        btmClip = 0.0f;
        int x1 = 0;
        int y1 = 0;
        int width = 0;
        int length = 0;
        int idx = 0;

        if(segmentsOn)
        {
            int segmentSize = heightPerInterval * dataHist.length / numSegments + 1;
            for(int i=0; i<numSegments; i++)
            {
                float colorNum = 1.0f - (float)i * 0.1f;
                g.setColor(new Color(colorNum,colorNum,colorNum));
                if(orientation == HORIZONTAL)
                {
                    g.fillRect(x + (i*segmentSize), y,  segmentSize, bwidth);
                }
                else
                {
                    g.fillRect(x, y + (i*segmentSize), bwidth, segmentSize);
                }
            }
        }

        if(verticalScale != -1.0f)
        {
            scale = bwidth / verticalScale;
        }

        for (int i = dataHist.length - 1; i >= 0; i--)
        {
            if (orientation == VERTICAL)
            {
                idx = 254 - i;
            }
            else
            {
                idx = i;

            }
            if(verticalScale != -1.0f)
            {
                if(dataHist[idx] > verticalScale)
                    dataHist[idx] = verticalScale;
            }
            if (idx > maxIdx)
            {
                g.setColor(Color.GRAY);
                topClip += dataHist[idx];
            }
            else if (idx < minIdx)
            {
                g.setColor(Color.GRAY);
                btmClip += dataHist[idx];
            }
            else
            {
                g.setColor(Color.BLUE);

                // draw a rectangle for each bar
            }
            if (orientation == VERTICAL)
            {
                x1 = x + bwidth - (int) (dataHist[idx] * scale);
                y1 = y + (i * heightPerInterval * (nPoints / dataHist.length));
                width = (int) (dataHist[idx] * scale);
                length = (heightPerInterval * (nPoints / dataHist.length));
            }
            else
            {
                x1 = x + (i * heightPerInterval * (nPoints / dataHist.length));
                y1 = y + bwidth - (int) (dataHist[idx] * scale);
                width = (heightPerInterval * (nPoints / dataHist.length));
                length = (int) (dataHist[idx] * scale);
            }
            g.fillRect(x1, y1, width, length);
        }
        g.setColor(Color.black);
        if (orientation == HORIZONTAL)
        {
            g.drawRect(x, y, heightPerInterval * dataHist.length, bwidth);
        }
        else
        {
            g.drawRect(x, y, bwidth, heightPerInterval * dataHist.length);

        }
        if (mouseBtn == 1)
        {
            paintOffScreen(offScreenGraphics);
            g.drawImage(offScreenImage, dataX, dataY /*- fm.getHeight()*/, this);
        }
        if (orientation != HORIZONTAL)
        {
            g.setColor(new Color(1.0f, 0.0f, 0.0f, 0.2f));
            g.fillPolygon(opacityX, opacityY, nPoints + 2);
            g.setColor(new Color(0.5f, 0.5f, 0.5f));
            g.drawPolyline(opacityX, opacityY, nPoints + 2);
        }
    }

    public void recalcClip()
    {
        int idx = 0;
        topClip = 0;
        btmClip = 0;
        for (int i = dataHist.length - 1; i >= 0; i--)
        {
            if (orientation == VERTICAL)
            {
                idx = 254 - i;
            }
            else
            {
                idx = i;

            }
            if (idx > maxIdx)
            {
                topClip += dataHist[idx];
            }
            else if (idx < minIdx)
            {
                btmClip += dataHist[idx];
            }
        }
    }

    private void addKeyPopup(int x, int y)
    {
        JPopupMenu tp = new JPopupMenu("Histogram Popup");
        this.add(tp);
        tp.add(kc);
        tp.show(this, x, y);
    }

    private void addMouseListener()
    {
        MouseListener mouseListener = new MouseAdapter()
        {
            public void mouseReleased(MouseEvent e)
            {
                dataValue = null;
                dataPercent = null;
                if (mouseBtn == 3)
                {
                    if (colorscale != null)
                    {
                        colorscalepanel.setOpacityValues(opacityValues);
                    }
                    mouseBtn = 0;
                    previousIdx = -1;
                    repaint();
                }
            }
            public void mousePressed(MouseEvent e)
            {
                mouseBtn = e.getButton();
                int mods = e.getModifiers();
                if((e.isShiftDown() && ((mods & InputEvent.BUTTON3_MASK)) != 0) || ((mods & InputEvent.BUTTON2_MASK) != 0))
                {
                    if(colorscale != null)
                        addKeyPopup(e.getX(),e.getY());
                }
            }
        };
        this.addMouseListener(mouseListener);
    }

    private void addMouseMotionListener()
    {
        MouseMotionListener motionListener = new MouseMotionAdapter()
        {
            public void mouseDragged(MouseEvent e)
            {
                int y = 0;
                int pickedY = e.getY();
                int pickedX = e.getX();

                if (orientation == HORIZONTAL)
                {
                    if (pickedY > barRect.y + barRect.height)
                        pickedY = barRect.y + barRect.height;
                    y = - (pickedX - xOrigin);
                }
                else
                {
                    if (pickedX > barRect.x + barRect.width)
                        pickedX = barRect.x + barRect.width;
                    y = pickedY - (barRect.y + barRect.height);

                }
                int index = -y / heightPerInterval;

                if (mouseBtn == e.BUTTON1)
                {
                    // Output the data value associated with the mouse position. SAJ
                    if ( (index < 0) || (index > maxIdx))
                    {
                        dataValue = null;
                        dataPercent = null;
                        return;
                    }
                    else
                    {
                        offScreenGraphics.dispose();
                        dataPercent = new Float(colorFormat.format(dataHist[index]));
                        dataValue = new Float(colorFormat.format(getDataValue(index)));

                        dataX = e.getX();
                        if ( (dataWidth + dataX) > barRect.width)
                        {
                            dataX = barRect.width - dataWidth;
                        }
                        else if (dataX < barRect.x)
                        {
                            dataX = barRect.x + 1;

                        }
                        dataY = e.getY() - fm.getHeight();
                        if ( (e.getY() - fm.getHeight()) < barRect.y)
                        {
                            dataY = barRect.y + 3;
                        }
                        else if (e.getY() > barRect.y + barRect.height)
                        {
                            dataY = barRect.y + barRect.height - fm.getHeight() - 3;
                        }
                    }
                    repaint();
                }
                else if ( (mouseBtn == e.BUTTON3) && (orientation != HORIZONTAL))
                {
                    if((index < minIdx) || (index > maxIdx))
                        return;

                    dataX = e.getX() + barRect.x;
                    if (dataX > barRect.width + barRect.x)
                        dataX = barRect.width + barRect.x;
                    else if (dataX < barRect.x)
                        dataX = barRect.x;

                    dataY = e.getY();
                    if (dataY < barRect.y)
                        dataY = barRect.y;
                    else if (dataY > barRect.y + barRect.height)
                        dataY = barRect.y + barRect.height;

                    opacityX[index + 1] = dataX;
                    opacityY[index + 1] = dataY;
                    opacityValues[index] = (float) (dataX - barRect.x) / barWidth;
                    if(opacityValues[index] > 1.0) opacityValues[index] = 1.0f;
                    if(opacityValues[index] < 0.0f) opacityValues[index] = 0.0f;

                        // Interpolated between previous index and current one
                    if((previousIdx != -1) && (Math.abs(previousIdx - index) > 1))
                    {
                        int length = previousIdx - index;
                        float scale = (float)((float)(opacityX[index+1] - opacityX[previousIdx + 1])/Math.abs(length));
                        if(length > 0)
                            for(int n = 1; n < length; n++)
                            {
                                opacityX[index + n + 1] = opacityX[index + 1] + (int)(scale * (float)n);
                                opacityY[index + n + 1] = opacityY[index + 1] + -(heightPerInterval * n);
                                opacityValues[index + n] = (float)(((float)opacityX[index + 1] + (scale * (float)n) - barRect.x) / barWidth);
                                if(opacityValues[index+n] > 1.0) opacityValues[index+n] = 1.0f;
                                if(opacityValues[index+n] < 0.0f) opacityValues[index+n] = 0.0f;
                            }
                        else
                            for(int n = -1; n > length; n--)
                            {
                                opacityX[index + n + 1] = opacityX[index + 1] + (int)(scale * (float)n);
                                opacityY[index + n + 1] = opacityY[index + 1] + -(heightPerInterval * n);
                                opacityValues[index + n] = (float)(((float)opacityX[index + 1] + (scale * (float)n) - barRect.x) / barWidth);
                                if(opacityValues[index+n] > 1.0) opacityValues[index+n] = 1.0f;
                                if(opacityValues[index+n] < 0.0f) opacityValues[index+n] = 0.0f;
                            }
                    }
                    previousIdx = index;
                    repaint();
                }
            }
        };
        this.addMouseMotionListener(motionListener);
    }

        public float getDataValue(int idx)
        {
            float scale = (maxVal - minVal) / nPoints;
            float value = idx * scale + minVal;
            return value;
        }

        public int getNumberOfSegments() { return numSegments; }
        public void setNumberOfSegment(int num)
        {
            numSegments = num;
        }
        public boolean getSegmentsOn() { return segmentsOn; }
        public void setSegmentsOn(boolean bool)
        {
            segmentsOn = bool;
        }

        private void paintOffScreen(Graphics og)
        {
            if(dataValue == null) return;
            og.clearRect(0, 0, dataWidth, dataHeight);
            og.setColor(Color.red);
            og.drawRect(0, 0, dataWidth - 1, dataHeight - 1);
            og.setColor(Color.black);
            og.setFont(floatFont);
            int x = (int) ( (float) dataWidth * 0.05f);
//            int y = (int) ((float) dataHeight * 0.10f);
            int y = (int) ( (float) dataHeight - (float) dataHeight * 0.10f);
            og.drawString(dataPercent.toString(), x, y);
//            og.drawString(dataValue.toString(), x, y);
//            og.drawString(dataPercent.toString(), x, y - fm.getHeight());
        }

        public boolean setClip(int min, int max)
        {
            boolean changed = false;
            if (minIdx != min)
            {
                minIdx = min;
                changed = true;
            }
            if (maxIdx != max)
            {
                maxIdx = max;
                changed = true;
            }
            return changed;
        }

        public void clearAll()
        {
            minIdx = 0;
            maxIdx = 254;
            dataMax = -1000.0f;
            dataMin = 1000.0f;
            topClip = 0.0f;
            btmClip = 0.0f;

            dataHist = null;
        }

        public int getNumberIndices()
        {
            return barHeight;
        }

        public float[] getOpacityValues()
        {
            return opacityValues;
        }

        public float getTopPercentageClipped()
        {
            return topClip;
        }
        public float getBottomPercentageClipped()
        {
            return btmClip;
        }
        public float[] getDataHistogram()
        {
            return dataHist;
        }
        public int getHeightPerInterval()
        {return heightPerInterval;
        }
        public void setHeightPerInterval(int heightPerInterval)
        {this.heightPerInterval = heightPerInterval;
        }
        public int getBarWidth()
        {return barWidth;
        }
        public void setBarWidth(int barWidth)
        {this.barWidth = barWidth;
        }
        public int getNPoints()
        {return nPoints;
        }
        public void setNPoints(int nPoints)
        {this.nPoints = nPoints;
        }

        static public void main(String[] args)
        {
            float[] dataHist =
                {
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.99f, 0.92f, 0.81f, 0.71f, 0.61f, 0.5f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                .9f, 0.75f, 0.23f, 0.04f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                0.0f
            };

            try
            {
                StsHistogramPanelOld histogramPanel = new StsHistogramPanelOld(HORIZONTAL, null);
                StsToolkit.createDialog(histogramPanel);
            }
            catch (Exception e)
            {
                e.printStackTrace();
            }
        }

    }
