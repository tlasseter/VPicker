package com.Sts.UI;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

import com.Sts.DBTypes.*;

import javax.swing.*;
import java.awt.*;
import java.awt.font.*;
import java.awt.geom.*;
import java.awt.image.*;

public class StsSuperGatherConfigPanel extends JPanel
{
    private StsPreStackLineSet volume = null;
    private BufferedImage bimg;
    private Ellipse2D ellipse = null;
    private Line2D line = null;
    Graphics2D g2d = null;
    Point[] gathers = null;
    AffineTransform at, aT;

    public StsSuperGatherConfigPanel()
    {
        setBackground(Color.white);
    }

    public void setPreStackSeismic(StsPreStackLineSet volume)
    {
        this.volume = volume;
    }
    /*
     * draws the transformed image, the String describing the current
     * transform, and the current transformation factors.
     */
    public void drawDemo(int w, int h, Graphics2D g2)
    {
        g2d = g2;
        Point2D point = null;
        Point2D pointT = new Point();
        int nPoints = 1, nRows = 1, nCols = 1;

        if(volume == null) return;
        g2d.setColor(Color.blue);

        Font font = g2d.getFont();
        FontRenderContext frc = g2d.getFontRenderContext();
        TextLayout tl = new TextLayout("Super Gather Layout", font, frc);
        tl.draw(g2d, (float) (w/2-tl.getBounds().getWidth()/2),(float)(tl.getAscent()+tl.getDescent()));

        float xPos = (float) (((3*w)/4)-tl.getBounds().getWidth()/2);
        float yPos = (float)(h-(tl.getAscent()+tl.getDescent()));

        if(volume instanceof StsPreStackLineSet3d)
            tl = new TextLayout("Inline Direction", font, frc);
        else
            tl = new TextLayout("Along Line Direction", font, frc);

        tl.draw(g2d, xPos, yPos);

        xPos = xPos + (float)tl.getBounds().getWidth();
        yPos = yPos - tl.getAscent()/2;
        g2d.drawLine((int)(xPos + 10), (int)yPos, (int)(xPos + w/10), (int)yPos);

        xPos = xPos  + (float)(w/10.0f);
        g2d.drawLine((int)xPos, (int)yPos, (int)(xPos - w/20), (int)(yPos + tl.getAscent()/4));
        g2d.drawLine((int)xPos, (int)yPos, (int)(xPos - w/20), (int)(yPos - tl.getAscent()/4));

        switch(volume.superGatherProperties.getGatherType())
        {
            case StsPreStackLineSetClass.SUPER_CROSS:
                nCols = volume.superGatherProperties.getNSuperGatherCols();
                nRows = volume.superGatherProperties.getNSuperGatherRows();
                nPoints = nCols + nRows;
                break;
            case StsPreStackLineSet3dClass.SUPER_RECT:
                nCols = volume.superGatherProperties.getNSuperGatherCols();
                nRows = volume.superGatherProperties.getNSuperGatherRows();
                nPoints = nCols * nRows;
                break;
            default:
                break;
        }
        gathers = new Point[nPoints];

        float dw = w * 0.6f;
        float dh = h * 0.6f;
        float dwi = dw / nCols;
        float dhi = dh / nRows;
        if(dwi > dhi)
            dwi = dhi;
        else
            dhi = dwi;

        if(volume instanceof StsPreStackLineSet3d)
        {
            if(volume.getXInc() > volume.getYInc())
                dhi = dhi * (volume.getYInc() / volume.getXInc());
            else
                dwi = dwi * (volume.getXInc() / volume.getYInc());
        }

        float sw = w * 0.2f;
        float sh = h * 0.15f;

        // Create Rotation Transform
        float centerX = sw + dw/2.0f;
        float centerY = h - sh - dh/2;
        aT = g2d.getTransform();
        at = AffineTransform.getTranslateInstance(centerX, centerY);
        //at.rotate(Math.toRadians(psClass.getCurrentSeismicVolume().getAngle()));
        g2d.transform(at);

        // Gathers
        g2d.setColor(Color.black);
        float startX = -4f;
        if(nCols > 1)
            startX = -(dwi*(nCols-1))/2.0f - 4f;
        float startY = -8f;
        if(nRows > 1)
            startY = (dhi*(nRows-1))/2.0f - 8f;

        int pointCount = 0;
        g2d.setColor(Color.blue);

        for (int j = 0; j < nCols; j++)
        {
            for (int i = 0; i < nRows; i++)
            {
                if((volume.superGatherProperties.getGatherType() == volume.lineSetClass.SUPER_CROSS) && ((i != nRows/2) && (j != nCols/2)))
                    continue;

                point = new Point((int)(startX + (dwi * j)),(int)(startY + (-dhi * i)));
                gathers[pointCount] = new Point();
                at.transform(point, gathers[pointCount]);
                pointCount++;
                ellipse = new Ellipse2D.Float(startX + (dwi * j), startY + (-dhi * i), 8, 8);
                if((j == nCols/2) && (i == nRows/2))
                {
                    g2d.setColor(Color.red);
                    g2d.fill(ellipse);
                }
                else
                {
                    g2d.setColor(Color.black);
                    g2d.draw(ellipse);
                }
            }
        g2d.transform(aT);
        }
    }

    public Graphics2D createGraphics2D(int w, int h)
    {
        Graphics2D g2 = null;
        if (bimg == null || bimg.getWidth() != w || bimg.getHeight() != h)
            bimg = (BufferedImage) createImage(w, h);

        g2 = bimg.createGraphics();
        g2.setBackground(getBackground());
        g2.clearRect(0, 0, w, h);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                            RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_RENDERING,
                            RenderingHints.VALUE_RENDER_QUALITY);
        return g2;
    }

    public void paint(Graphics g)
    {
        Dimension d = getSize();
        Graphics2D g2 = createGraphics2D(d.width, d.height);
        drawDemo(d.width, d.height, g2);
        g2.dispose();
        g.drawImage(bimg, 0, 0, this);
    }

}
