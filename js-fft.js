/*
The MIT License

Copyright (c) 2012 Nick Jones

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

// Math constants and functions we need.
const PI = Math.PI;
const SQRT1_2 = Math.SQRT1_2;

class baseComplexArray {
  constructor(other, arrayType = Float32Array) {
    if (other instanceof ComplexArray) {
      // Copy constuctor.
      this.ArrayType = other.ArrayType;
      this.real = new this.ArrayType(other.real);
      this.imag = new this.ArrayType(other.imag);
    } else {
      this.ArrayType = arrayType;
      // other can be either an array or a number.
      this.real = new this.ArrayType(other);
      this.imag = new this.ArrayType(this.real.length);
    }

    this.length = this.real.length;
  }

  toString() {
    const components = [];

    this.forEach((value, i) => {
      components.push(
        `(${value.real.toFixed(2)}, ${value.imag.toFixed(2)})`
      );
    });

    return `[${components.join(', ')}]`;
  }

  forEach(iterator) {
    const n = this.length;
    // For gc efficiency, re-use a single object in the iterator.
    const value = Object.seal(Object.defineProperties({}, {
      real: { writable: true }, imag: { writable: true },
    }));

    for (let i = 0; i < n; i++) {
      value.real = this.real[i];
      value.imag = this.imag[i];
      iterator(value, i, n);
    }
  }

  // In-place mapper.
  map(mapper) {
    this.forEach((value, i, n) => {
      mapper(value, i, n);
      this.real[i] = value.real;
      this.imag[i] = value.imag;
    });

    return this;
  }

  conjugate() {
    return new ComplexArray(this).map((value) => {
      value.imag *= -1;
    });
  }

  toMagnitudeAndPhase() {
    const mags = [];

    this.forEach((value, i) => {
      mags.push(
        { magnitude: Math.sqrt(value.real * value.real + value.imag * value.imag), phase: Math.atan(value.imag, value.real) * 180 / PI }
      );
    })

    return mags;
  }

  toRealAndImag() {
    const components = [];

    this.forEach((value, i) => {
      components.push(
        { real: value.real, imag: value.imag }
      );
    })

    return components;
  }
}

class ComplexArray extends baseComplexArray {
  FFT() {
    return this.fft(this, false);
  }

  InvFFT() {
    return this.fft(this, true);
  }

  // Applies a frequency-space filter to input, and returns the real-space
  // filtered input.
  // filterer accepts freq, i, n and modifies freq.real and freq.imag.
  frequencyMap(filterer) {
    return this.FFT().map(filterer).InvFFT();
  }

  fft(input, inverse) {
    const n = input.length;

    if (n & (n - 1)) {
      return this.FFT_Recursive(input, inverse);
    } else {
      return this.FFT_2_Iterative(input, inverse);
    }
  }

  FFT_Recursive(input, inverse) {
    const n = input.length;

    if (n === 1) {
      return input;
    }

    const output = new ComplexArray(n, input.ArrayType);

    // Use the lowest odd factor, so we are able to use FFT_2_Iterative in the
    // recursive transforms optimally.
    const p = this.LowestOddFactor(n);
    const m = n / p;
    const normalisation = 1 / Math.sqrt(p);
    let recursive_result = new ComplexArray(m, input.ArrayType);

    // Loops go like O(n Σ p_i), where p_i are the prime factors of n.
    // for a power of a prime, p, this reduces to O(n p log_p n)
    for (let j = 0; j < p; j++) {
      for (let i = 0; i < m; i++) {
        recursive_result.real[i] = input.real[i * p + j];
        recursive_result.imag[i] = input.imag[i * p + j];
      }
      // Don't go deeper unless necessary to save allocs.
      if (m > 1) {
        recursive_result = this.fft(recursive_result, inverse);
      }

      const del_f_r = Math.cos(2 * PI * j / n);
      const del_f_i = (inverse ? -1 : 1) * Math.sin(2 * PI * j / n);
      let f_r = 1;
      let f_i = 0;

      for (let i = 0; i < n; i++) {
        const _real = recursive_result.real[i % m];
        const _imag = recursive_result.imag[i % m];

        output.real[i] += f_r * _real - f_i * _imag;
        output.imag[i] += f_r * _imag + f_i * _real;

        [f_r, f_i] = [
          f_r * del_f_r - f_i * del_f_i,
          f_i = f_r * del_f_i + f_i * del_f_r,
        ];
      }
    }

    // Copy back to input to match FFT_2_Iterative in-placeness
    // TODO: faster way of making this in-place?
    for (let i = 0; i < n; i++) {
      input.real[i] = normalisation * output.real[i];
      input.imag[i] = normalisation * output.imag[i];
    }

    return input;
  }

  FFT_2_Iterative(input, inverse) {
    const n = input.length;

    const output = this.BitReverseComplexArray(input);
    const output_r = output.real;
    const output_i = output.imag;
    // Loops go like O(n log n):
    //   width ~ log n; i,j ~ n
    let width = 1;
    while (width < n) {
      const del_f_r = Math.cos(PI / width);
      const del_f_i = (inverse ? -1 : 1) * Math.sin(PI / width);
      for (let i = 0; i < n / (2 * width); i++) {
        let f_r = 1;
        let f_i = 0;
        for (let j = 0; j < width; j++) {
          const l_index = 2 * i * width + j;
          const r_index = l_index + width;

          const left_r = output_r[l_index];
          const left_i = output_i[l_index];
          const right_r = f_r * output_r[r_index] - f_i * output_i[r_index];
          const right_i = f_i * output_r[r_index] + f_r * output_i[r_index];

          output_r[l_index] = SQRT1_2 * (left_r + right_r);
          output_i[l_index] = SQRT1_2 * (left_i + right_i);
          output_r[r_index] = SQRT1_2 * (left_r - right_r);
          output_i[r_index] = SQRT1_2 * (left_i - right_i);

          [f_r, f_i] = [
            f_r * del_f_r - f_i * del_f_i,
            f_r * del_f_i + f_i * del_f_r,
          ];
        }
      }
      width <<= 1;
    }

    return output;
  }

  BitReverseIndex(index, n) {
    let bitreversed_index = 0;

    while (n > 1) {
      bitreversed_index <<= 1;
      bitreversed_index += index & 1;
      index >>= 1;
      n >>= 1;
    }
    return bitreversed_index;
  }

  BitReverseComplexArray(array) {
    const n = array.length;
    const flips = new Set();

    for (let i = 0; i < n; i++) {
      const r_i = this.BitReverseIndex(i, n);

      if (flips.has(i)) continue;

      [array.real[i], array.real[r_i]] = [array.real[r_i], array.real[i]];
      [array.imag[i], array.imag[r_i]] = [array.imag[r_i], array.imag[i]];

      flips.add(r_i);
    }

    return array;
  }

  LowestOddFactor(n) {
    const sqrt_n = Math.sqrt(n);
    let factor = 3;

    while (factor <= sqrt_n) {
      if (n % factor === 0) return factor;
      factor += 2;
    }
    return n;
  }

}

module.exports = function (RED) {
  function JsFftNode(config) {

    RED.nodes.createNode(this, config);
    var node = this;
    node.algorithm = config.algorithm || "fft";
    node.freqmin = config.freqmin;
    node.freqmax = config.freqmax;

    node.on('input', function (msg) {
      var data = new ComplexArray(msg.payload.length);
      switch (node.algorithm) {
        case "fft":
          data.map((value, i, n) => {
            value.real = msg.payload[i];
          });
          data.FFT();
          msg.payload = data.toRealAndImag();
          break;s

        case "inv":
          data.map((value, i, n) => {
            value.real = msg.payload[i].real;
            value.imag = msg.payload[i].imag;
          });
          data.InvFFT();
          msg.payload = data.toRealAndImag();
          break;

        case "mgn":
          data.map((value, i, n) => {
            value.real = msg.payload[i];
          });
          data.FFT();
          msg.payload = data.toMagnitudeAndPhase();
          break;

        default:
          msg.payload = "Type not specified";
          break;
      }

      node.send(msg);
    });
  }
  RED.nodes.registerType("js-fft", JsFftNode);
}