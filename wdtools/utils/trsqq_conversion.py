from random import sample
import re


class TrsqqConverter:
    def __init__(self):
        self.directions = list('NSEW')

    def convert(self, x):
        x = f'{x:<08s}'
        xt = (
            self._get_tr_code(x)
            + self._get_tr_code(x)
            + self._get_s_code(x)
            + self._get_qq_code(x))
        return xt[:16]

    def _get_tr_code(self, x):
        'Get township/range code from the taxlot'
        nms = re.findall('\d+', x)
        nm1 = nms[1][0:2]
        lts = re.findall("[a-zA-Z]+", x)
        if lts[1] in self.directions:
            lts2 = lts[1]
        else:
            lts2 = sample(self.directions, 1)[0]
        if len(nm1) == 1:
            tr1 = f'0{nm1}'
        elif len(nm1) == 3:
            tr1 = nm1[1:3]
        else:
            tr1 = nm1
        if ('V' in lts2) or ('Y' in lts2):
            tr2 = '.50'
            tr3 = lts2[1]
        elif ('X' in lts2) or ('Z' in lts2):
            if any([x in lts2 for x in ['XS', 'ZN', 'XE', 'ZW']]):
                tr2 = '.75'
            else:
                tr2 = '.25'
        else:
            tr2 = '.00'
            tr3 = lts2
        res = tr1 + tr2 + tr3
        return res

    @staticmethod
    def _get_s_code(x):
        'Get section code from trsqq code'
        if len(x) <= 7:
            s = '00'
        else:
            nms = re.findall('\d+', x)
            n = len(nms)
            k = len(nms[1])
            if (n < 3) and (k > 2):
                nm1 = nms[1][(k-2):(k+1)]
            else:
                nm1 = nms[2]
            if len(nm1) == 1:
                s = f'0{nm1}'
            else:
                s = nm1
        return s

    def _get_qq_code(self, x):
        'Get QQ code from trsqq code'
        nms = re.findall('\d+', x)
        lts = re.findall("[a-zA-Z]+", x)
        if (len(lts) == 2) and (lts[1] not in self.directions):
            if len(lts[1]) == 2:
                qq = lts[1]
            else:
                qq = f'{lts[1]}0'
        elif len(nms[2]) > 2:
            if nms[2] == 4:
                qq = nms[2][2:4]
            else:
                qq = f'{nms[2][2]}0'
        elif len(lts) == 3:
            if len(lts[2]) == 2:
                qq = lts[2]
            else:
                qq = f'{lts[2]}0'
        else:
            if len(nms[0]) == 1:
                t = f'0{nms[0]}'
            else:
                t = nms[0]
            if len(nms[1]) == 1:
                r = f'0{nms[1]}'
            else:
                r = nms[1]
            if len(nms[2]) == 1:
                s = f'0{nms[2]}'
            else:
                s = nms[2]
            trsqq = t + lts[0] + r + lts[1] + s
            qq = f'{trsqq:0<10}'[8:10]
        return qq
