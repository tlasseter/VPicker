package com.Sts.UI.Beans;

/**
 * ArrowRoller.java  1.00 97/07/09 Merlin Hughes
 * <p/>
 * Copyright (c) 1997 Merlin Hughes, All Rights Reserved.
 * <p/>
 * Permission to use, copy, modify, and distribute this software
 * for commercial and non-commercial purposes and without fee is
 * hereby granted provided that this copyright notice appears in
 * all copies.
 * <p/>
 * http://prominence.com/                         ego@merlin.org
 */

class ArrowRoller extends Thread
{
    private static final int DELAY_MS = 400, REPEAT_MS = 100;

    ArrowRoller()
    {
        super("ArrowRoller");
        start();
    }

    private ArrowBean target;
    private int state;

    synchronized void abort()
    {
        stop();
        notify();
        target = null;
        ++state;
    }

    synchronized void addTarget(ArrowBean t)
    {
        target = t;
        ++state;
        notify();
    }

    synchronized void removeTarget(ArrowBean t)
    {
        if (target == t)
        {
            target = null;
            ++state;
        }
    }

    public void run()
    {
        while (true)
        {
            try
            {
                ArrowBean target;
                int state;
                synchronized (this)
                {
                    while ((target = this.target) == null)
                        wait();
                    state = this.state;
                }
                sleep(DELAY_MS);
                while (state == this.state)
                {
                    target.fireActionEvent();
                    sleep(REPEAT_MS);
                }
            }
            catch (InterruptedException ex)
            {
            }
        }
    }
}
