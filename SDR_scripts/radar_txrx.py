#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# Author: nuc1
# GNU Radio version: v3.8.5.0-5-g982205bd

from datetime import datetime
from gnuradio import blocks
import pmt
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time


class radar_txrx(gr.top_block):

    def __init__(self, duration=30, freq=1.65e9, input_file='/home/nuc1/uav-bistatic-sar/SDR_scripts/tx_wave_56.dat', rx_gain=40, samp_rate=40e6, tx_gain=40):
        gr.top_block.__init__(self, "Not titled yet")

        ##################################################
        # Parameters
        ##################################################
        self.duration = duration
        self.freq = freq
        self.input_file = input_file
        self.rx_gain = rx_gain
        self.samp_rate = samp_rate
        self.tx_gain = tx_gain

        ##################################################
        # Variables
        ##################################################
        self.output_dir = output_dir = "/home/nuc1/uav-bistatic-sar/SDR_scripts/"
        self.output_file = output_file = output_dir+datetime.now().strftime("%Y-%m-%d %H-%M-%S")+ "_f"+ str(int(freq / 1e6))+ "_s" + str(int(samp_rate / 1e6))+ ".dat"

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("", "num_recv_frames=1024")),
            uhd.stream_args(
                cpu_format="sc16",
                args='',
                channels=list(range(0,1)),
            ),
        )
        self.uhd_usrp_source_0.set_time_source('gpsdo', 0)
        self.uhd_usrp_source_0.set_clock_source('gpsdo', 0)
        self.uhd_usrp_source_0.set_center_freq(freq, 0)
        self.uhd_usrp_source_0.set_rx_agc(False, 0)
        self.uhd_usrp_source_0.set_gain(rx_gain, 0)
        self.uhd_usrp_source_0.set_antenna('RX2', 0)
        self.uhd_usrp_source_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_0.set_min_output_buffer(100000000)
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
            ",".join(("", "")),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
            '',
        )
        self.uhd_usrp_sink_0.set_center_freq(freq, 0)
        self.uhd_usrp_sink_0.set_gain(tx_gain, 0)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_sink_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_sink_0.set_min_output_buffer(100000000)
        self.blocks_stream_to_vector_0 = blocks.stream_to_vector(gr.sizeof_short*2, 1)
        self.blocks_head_0_0 = blocks.head(gr.sizeof_gr_complex*1, int(samp_rate *duration*1.01))
        self.blocks_head_0 = blocks.head(gr.sizeof_short*2, int(samp_rate *duration))
        self.blocks_file_source_0_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, input_file, True, 0, 0)
        self.blocks_file_source_0_0_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_sink_1 = blocks.file_sink(gr.sizeof_short*2, output_file, False)
        self.blocks_file_sink_1.set_unbuffered(False)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_0_0_0, 0), (self.blocks_head_0_0, 0))
        self.connect((self.blocks_head_0, 0), (self.blocks_file_sink_1, 0))
        self.connect((self.blocks_head_0_0, 0), (self.uhd_usrp_sink_0, 0))
        self.connect((self.blocks_stream_to_vector_0, 0), (self.blocks_head_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.blocks_stream_to_vector_0, 0))


    def get_duration(self):
        return self.duration

    def set_duration(self, duration):
        self.duration = duration
        self.blocks_head_0.set_length(int(self.samp_rate *self.duration))
        self.blocks_head_0_0.set_length(int(self.samp_rate *self.duration*1.01))

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq
        self.set_output_file(self.output_dir+datetime.now().strftime("%Y-%m-%d %H-%M-%S")+ "_f"+ str(int(self.freq / 1e6))+ "_s" + str(int(self.samp_rate / 1e6))+ ".dat")
        self.uhd_usrp_sink_0.set_center_freq(self.freq, 0)
        self.uhd_usrp_source_0.set_center_freq(self.freq, 0)

    def get_input_file(self):
        return self.input_file

    def set_input_file(self, input_file):
        self.input_file = input_file
        self.blocks_file_source_0_0_0.open(self.input_file, True)

    def get_rx_gain(self):
        return self.rx_gain

    def set_rx_gain(self, rx_gain):
        self.rx_gain = rx_gain
        self.uhd_usrp_source_0.set_gain(self.rx_gain, 0)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_output_file(self.output_dir+datetime.now().strftime("%Y-%m-%d %H-%M-%S")+ "_f"+ str(int(self.freq / 1e6))+ "_s" + str(int(self.samp_rate / 1e6))+ ".dat")
        self.blocks_head_0.set_length(int(self.samp_rate *self.duration))
        self.blocks_head_0_0.set_length(int(self.samp_rate *self.duration*1.01))
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 0)
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_bandwidth(self.samp_rate, 0)

    def get_tx_gain(self):
        return self.tx_gain

    def set_tx_gain(self, tx_gain):
        self.tx_gain = tx_gain
        self.uhd_usrp_sink_0.set_gain(self.tx_gain, 0)

    def get_output_dir(self):
        return self.output_dir

    def set_output_dir(self, output_dir):
        self.output_dir = output_dir
        self.set_output_file(self.output_dir+datetime.now().strftime("%Y-%m-%d %H-%M-%S")+ "_f"+ str(int(self.freq / 1e6))+ "_s" + str(int(self.samp_rate / 1e6))+ ".dat")

    def get_output_file(self):
        return self.output_file

    def set_output_file(self, output_file):
        self.output_file = output_file
        self.blocks_file_sink_1.open(self.output_file)




def argument_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-d", "--duration", dest="duration", type=eng_float, default="30.0",
        help="Set duration in seconds [default=%(default)r]")
    parser.add_argument(
        "-f", "--freq", dest="freq", type=eng_float, default="1.65G",
        help="Set freq [default=%(default)r]")
    parser.add_argument(
        "-i", "--input-file", dest="input_file", type=str, default='/home/nuc1/uav-bistatic-sar/SDR_scripts/tx_wave_56.dat',
        help="Set input_file [default=%(default)r]")
    parser.add_argument(
        "--rx-gain", dest="rx_gain", type=eng_float, default="40.0",
        help="Set rx_gain [default=%(default)r]")
    parser.add_argument(
        "-r", "--samp-rate", dest="samp_rate", type=eng_float, default="40.0M",
        help="Set samp_rate [default=%(default)r]")
    parser.add_argument(
        "--tx-gain", dest="tx_gain", type=eng_float, default="40.0",
        help="Set tx_gain [default=%(default)r]")
    return parser


def main(top_block_cls=radar_txrx, options=None):
    if options is None:
        options = argument_parser().parse_args()
    tb = top_block_cls(duration=options.duration, freq=options.freq, input_file=options.input_file, rx_gain=options.rx_gain, samp_rate=options.samp_rate, tx_gain=options.tx_gain)

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    tb.wait()


if __name__ == '__main__':
    main()
