package com.Sts.UI;


/**
 * <p>Title: </p>
 * <p>Description: Extension of FlowLayout which correctly resizes</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

/** FlowLayout doesn't properly resize.  So we have overridden the layoutContainer
 *  method, setting the container height to the minimum size needed to accommodate
 *  all the components.
 */

import java.awt.*;

public class StsFlowLayout extends FlowLayout
{

    public StsFlowLayout()
    {
    }

    public StsFlowLayout(int align)
    {
	    this(align, 5, 5);
    }

    public StsFlowLayout(int align, int hgap, int vgap)
    {
        super(align, hgap, vgap);
    }


    public void layoutContainer(Container target)
    {
        synchronized (target.getTreeLock())
        {
            int width = target.getWidth();
            if(width == 0) return;
            Dimension size = getOptimumSize(target);
            target.setSize(size);
        }
    }

    /** Using current target width, layout member components in rows, wrapping
     *  to next row if width is exceeded.
     */
    private Dimension getOptimumSize(Container target)
    {
        int hgap = getHgap();
        int vgap = getVgap();

        Insets insets = target.getInsets();
        int maxwidth = target.getWidth() - (insets.left + insets.right + hgap*2);
        int nmembers = target.getComponentCount();
        int x = 0, y = insets.top + vgap;
        int rowh = 0, start = 0;

        boolean ltr = target.getComponentOrientation().isLeftToRight();

        for (int i = 0 ; i < nmembers ; i++)
        {
            Component m = target.getComponent(i);
            if (m.isVisible())
            {
                Dimension d = m.getPreferredSize();
                m.setSize(d.width, d.height);

                if ((x == 0) || ((x + d.width) <= maxwidth))
                {
                    if (x > 0)
                    {
                        x += hgap;
                    }
                    x += d.width;
                    rowh = Math.max(rowh, d.height);
                }
                else
                {
                    moveComponents(target, insets.left + hgap, y, maxwidth - x, rowh, start, i, ltr);
                    x = d.width;
                    y += vgap + rowh;
                    rowh = d.height;
                    start = i;
                }
            }
        }
        moveComponents(target, insets.left + hgap, y, maxwidth - x, rowh, start, nmembers, ltr);

//        Dimension size = target.getSize();
//        System.out.println("size: " + size.toString());

        Dimension newSize = new Dimension(target.getWidth(), y+rowh+vgap+insets.bottom);
//        System.out.println("newSize: " + newSize.toString());
        return newSize;
    }

    private void moveComponents(Container target, int x, int y, int width, int height,
                                int rowStart, int rowEnd, boolean ltr)
    {
        synchronized (target.getTreeLock())
        {
            switch (getAlignment())
            {
                case LEFT:
                    x += ltr ? 0 : width;
                    break;
                case CENTER:
                    x += width / 2;
                    break;
                case RIGHT:
                    x += ltr ? width : 0;
                    break;
                case LEADING:
                    break;
                case TRAILING:
                    x += width;
                    break;
            }
            for (int i = rowStart ; i < rowEnd ; i++)
            {
                Component m = target.getComponent(i);

                if (m.isVisible())
                {
                    int mWidth = m.getWidth();
                    int mHeight = m.getHeight();
                    if (ltr)
                    {
                        m.setLocation(x, y + (height - mHeight) / 2);
                    }
                    else
                    {
                        m.setLocation(target.getWidth() - x - mWidth, y + (height - mHeight) / 2);
                    }
                    x += mWidth + getHgap();
                }
            }
        }
    }

    public Dimension preferredLayoutSize(Container target)
    {
        int width = target.getWidth();
        if(width == 0) return super.preferredLayoutSize(target);

        synchronized (target.getTreeLock())
        {
            return getOptimumSize(target);
        }
    }
}