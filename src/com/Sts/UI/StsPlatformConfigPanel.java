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
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.*;
import java.awt.geom.*;
import java.awt.image.*;

public class StsPlatformConfigPanel extends JPanel implements MouseListener
{
    private BufferedImage bimg;
    private StsPlatform platform;
    private Ellipse2D ellipse = null;
    private Line2D line = null;
    private boolean labelSlots = false;
    Graphics2D g2d = null;
    Point[] slots = null;
    int slotIdx = -1;
    AffineTransform at, aT;

    public StsPlatformConfigPanel()
    {
        setBackground(Color.white);
        addMouseListener(this);
    }

    public void setPlatform(StsPlatform platform)
    {
        this.platform = platform;
    }

    public void setLabelSlots(boolean val)
    {
        labelSlots = val;
    }
    public boolean getLabelSlots() { return labelSlots; }

    /*
     * draws the transformed image, the String describing the current
     * transform, and the current transformation factors.
     */
    public void drawDemo(int w, int h, Graphics2D g2)
    {
        g2d = g2;
        Point2D point = null;
        Point2D pointT = new Point();

        if(platform == null)
            return;
        g2d.setColor(Color.blue);
        Font font = g2d.getFont();
        FontRenderContext frc = g2d.getFontRenderContext();
        TextLayout tl = new TextLayout(platform.getName(), font, frc);
        tl.draw(g2d, (float) (w/2-tl.getBounds().getWidth()/2),(float) (tl.getAscent()+tl.getDescent()));
        slots = new Point[platform.getNumSlots()];

        if(platform.checkValidity())
        {
            if(platform.getSlotLayoutByString().equals(platform.SLOT_LAYOUT_STRINGS[platform.RECTANGULAR]))
            {
                float dw = w * 0.6f;
                float dh = h * 0.6f;
                float dwi = dw / platform.getNCols();
                float dhi = dh / platform.getNRows();
                if(dwi > dhi)
                    dwi = dhi;
                else
                    dhi = dwi;
                if(platform.getColSpacing() > platform.getRowSpacing())
                    dhi = dhi * (platform.getRowSpacing() / platform.getColSpacing());
                else
                    dwi = dwi * (platform.getColSpacing() / platform.getRowSpacing());
                float sw = w * 0.2f;
                float sh = h * 0.15f;

                // Origin Text
                g2d.setColor(Color.blue);
                String originStg = new String(String.valueOf("(" + platform.getXOrigin() + ", " + String.valueOf(platform.getYOrigin()) + ")"));
                tl = new TextLayout(originStg, font, frc);
                tl.draw(g2d, (float)(w/2-(tl.getBounds().getWidth()/2.0f)), h + (tl.getAscent()+tl.getDescent()*2));

                // Vertical North line
                g2d.setColor(Color.red);
                line = new Line2D.Float(w *0.05f, h*0.10f, w *0.05f, h*0.05f );
                g2d.draw(line);
                line = new Line2D.Float(w *0.06f, h*0.10f, w *0.06f, h*0.05f );
                g2d.draw(line);
                line = new Line2D.Float(w *0.045f, h*0.07f, w *0.055f, h*0.03f );
                g2d.draw(line);
                line = new Line2D.Float(w *0.065f, h*0.07f, w *0.055f, h*0.03f );
                g2d.draw(line);
                tl = new TextLayout("North", font, frc);
                tl.draw(g2d, (float)(w*0.055f -tl.getBounds().getWidth()/2.0f) , (float)(h*.10f + (tl.getAscent()+tl.getDescent()*2)));

                // Create Rotation Transform
                float centerX = sw + dw/2.0f;
                float centerY = h - sh - dh/2;
                aT = g2d.getTransform();
                at = AffineTransform.getTranslateInstance(centerX, centerY);
                at.rotate(Math.toRadians(platform.getRotationAngle()));
                g2d.transform(at);

                // Draw center of pad
                g2d.setColor(Color.blue);
                line = new Line2D.Float(2, 0, 5, 0);
                g2d.draw(line);
                line = new Line2D.Float(-2, 0, -5, 0);
                g2d.draw(line);
                line = new Line2D.Float(0, 2, 0, 5 );
                g2d.draw(line);
                line = new Line2D.Float(0, -2, 0, -5 );
                g2d.draw(line);

                // Rotating angle line
                line = new Line2D.Float(0, -dh/2, 0, -(dh/2 + h*0.075f));
                g2d.draw(line);
                line = new Line2D.Float(w*0.01f, -dh/2, w*0.01f, -(dh/2 + h*0.075f));
                g2d.draw(line);
                line = new Line2D.Float(-w*0.01f, -(dh/2+h*0.055f), w*0.005f, -(dh/2 + h*0.10f));
                g2d.draw(line);
                line = new Line2D.Float(w*0.02f, -(dh/2+h*0.055f), w*0.005f, -(dh/2 + h*0.10f));
                g2d.draw(line);

                // Rotation Angle Text
                g2d.setColor(Color.blue);
                tl = new TextLayout(String.valueOf(platform.getRotationAngle()), font, frc);
                tl.draw(g2d, (float) -tl.getBounds().getWidth()/2.0f, -(dh/2 + h*0.10f)-tl.getDescent());

                // Slots
                g2d.setColor(Color.black);
                float startX = -4f;
                if(platform.getNCols() > 1)
                    startX = -(dwi*(platform.getNCols()-1))/2.0f - 4f;
                float startY = -4f;
                if(platform.getNRows() > 1)
                    startY = (dhi*(platform.getNRows()-1))/2.0f - 4f;

                for (int j = 0; j < platform.getNRows(); j++)
                    for (int i = 0; i < platform.getNCols(); i++)
                    {
                        point = new Point((int)(startX + (dwi * i)),(int)(startY + (-dhi * j)));
                        slots[j*platform.getNCols() + i] = new Point();
                        at.transform(point, slots[j*platform.getNCols() + i]);

                        ellipse = new Ellipse2D.Float(startX + (dwi * i), startY + (-dhi * j), 8, 8);
                        if(platform.getCurrentSlotIndex() == j*platform.getNCols() + i)
                            g2d.setColor(Color.red);
                        else
                            g2d.setColor(Color.black);
                        g2d.draw(ellipse);
                        if(platform.isWellInSlot(j*platform.getNCols() + i))
                        {
                            // Draw well trajectory arrows
                            g2d.fill(ellipse);
                            drawWellArrow(point, platform.getSlotName(j*platform.getNCols() + i));
                        }

                        if(labelSlots)
                        {
                            tl = new TextLayout(platform.getSlotName(j*platform.getNCols() + i), font, frc);
                            tl.draw(g2d, (float) startX + (dwi * i), startY + (-dhi * j));
                        }
                    }
                g2d.transform(aT);
            }
            else
            {
                // Circular Layout
                g2d.setColor(Color.red);
                tl = new TextLayout("Circular Layout is not currently supportted", font, frc);
                tl.draw(g2d, (float)(w/2-tl.getBounds().getWidth()/2),(float)(h/2));
            }
        }
        else
        {
            TextLayout t2 = new TextLayout("Incomplete Definition", font, frc);
            g2d.setColor(Color.red);
            t2.draw(g2d, (float) (w / 2 - t2.getBounds().getWidth() / 2), (float) (h /2));
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

    public void mouseExited(MouseEvent m) {}
    public void mouseEntered(MouseEvent m) {}
    public void mouseReleased(MouseEvent m)
    {
        Point2D point = m.getPoint();

        for(int i=0; i< slots.length; i++)
        {
            if((Math.abs(point.getX()-slots[i].getX()) < 8) && (Math.abs(point.getY()-slots[i].getY()) < 8))
            {
                platform.setCurrentSlotIndex(i);
                repaint();
                break;
            }
        }
    }
    public void mouseClicked(MouseEvent m) {}
    public void mousePressed(MouseEvent m)
    {
        platform.setCurrentSlotIndex(-1);
    }

    public void drawWellArrow(Point2D slotCenter, String slotName)
    {
        float angleD = platform.getWellDirectionAngle(slotName);
        if(angleD == StsParameters.nullValue) return;
 //       StsWell well = platform.getWellAtSlot(slotName);
//        double relativeX = well.getXOrigin() - well.getCurrentModel().getProject().getXOrigin();
//        double relativeY = well.getYOrigin() - well.getCurrentModel().getProject().getYOrigin();

        // Calculate top hole to bottom hole angle

        double delX = 20*StsMath.sind(angleD - platform.getRotationAngle());
        double delY = 20*StsMath.cosd(angleD - platform.getRotationAngle());
   /*
        double deltaX = (double)(relativeX - well.getPointAtMDepth(well.getMaxMDepth(), false).getX());
        double deltaY = (double)(relativeY - well.getPointAtMDepth(well.getMaxMDepth(), false).getY());
        double delta = Math.sqrt((deltaX * deltaX) + (deltaY * deltaY));
        double angle = Math.atan(deltaX/deltaY);
        double delX = Math.sin(angle- StsMath.RADperDEG*platform.getRotationAngle())*20;
        double delY = Math.cos(angle- StsMath.RADperDEG*platform.getRotationAngle())*20;
    */
        Line2D line = new Line2D.Float((float)slotCenter.getX() + 4, (float)slotCenter.getY() + 4,
                                       (float)(slotCenter.getX() + 4 + delX), (float)(slotCenter.getY() + 4 - delY));
        g2d.draw(line);
    }
/*
    public void drawWellPath(Point2D slotCenter, String slotName)
    {
        double deltaX, deltaY, angle, delX, delY;

        StsWell well = platform.getWellAtSlot(slotName);
        double relativeX = well.getXOrigin() - well.getCurrentModel().getProject().getXOrigin();
        double relativeY = well.getYOrigin() - well.getCurrentModel().getProject().getYOrigin();

        Line2D line = new Line2D.Float();

        float slotXOrigin = (float)slotCenter.getX() + 4.0f;
        float slotYOrigin = (float)slotCenter.getY() + 4.0f;
        float ex, ey, sx, sy;
        sx = slotXOrigin;
        sy = slotYOrigin;
        StsPoint[] points = well.getPoints();
//        float xScale = ((float)getWidth()/2.0f)/(float)points.length;
//        float yScale = ((float)getHeight()/2.0f)/(float)points.length;

        for(int i=0; i<points.length; i=i+5)
        {
            angle = 0.0;
            deltaX = (double)(relativeX - points[i].getX());
            deltaY = (double)(relativeY - points[i].getY());
            if(deltaY != 0.0)
                angle = Math.atan(deltaX/deltaY);

            delX = Math.sin(angle- StsMath.RADperDEG*platform.getRotationAngle())*5.0f;
            delY = Math.cos(angle- StsMath.RADperDEG*platform.getRotationAngle())*5.0f;

            ex = (float)(sx + delX);
            ey = (float)(sy - delY);

            line.setLine(sx, sy, ex, ey);
            g2d.draw(line);

            sx = ex;
            sy = ey;
        }
    }
*/

}
