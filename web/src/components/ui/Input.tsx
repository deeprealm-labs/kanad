/**
 * Professional Input Component - IBM-Inspired Design
 *
 * Features:
 * - Sharp rectangular design (no rounded corners)
 * - IBM-inspired color palette
 * - Multiple sizes (sm, md, lg)
 * - State variants (error, success)
 * - Accessible with labels and error messages
 */

import React from 'react';

export interface InputProps extends React.InputHTMLAttributes<HTMLInputElement> {
  label?: string;
  error?: string;
  success?: boolean;
  helpText?: string;
  inputSize?: 'sm' | 'md' | 'lg';
  fullWidth?: boolean;
}

const Input = React.forwardRef<HTMLInputElement, InputProps>(
  (
    {
      label,
      error,
      success,
      helpText,
      inputSize = 'md',
      fullWidth = true,
      className = '',
      id,
      ...props
    },
    ref
  ) => {
    // Generate unique ID if not provided
    const inputId = id || `input-${Math.random().toString(36).substr(2, 9)}`;

    // Base classes
    const baseClasses = 'sharp transition-fast';

    // Size classes
    const sizeClasses = {
      sm: 'input-sm',
      md: '', // Default from global styles
      lg: 'input-lg',
    };

    // State classes
    const stateClasses = error
      ? 'input-error'
      : success
      ? 'input-success'
      : '';

    // Width classes
    const widthClasses = fullWidth ? 'w-full' : '';

    // Combine all classes
    const inputClasses = `
      ${baseClasses}
      ${sizeClasses[inputSize]}
      ${stateClasses}
      ${widthClasses}
      ${className}
    `.trim().replace(/\s+/g, ' ');

    return (
      <div className={fullWidth ? 'w-full' : ''}>
        {label && (
          <label htmlFor={inputId} className="block mb-2">
            {label}
          </label>
        )}
        <input
          ref={ref}
          id={inputId}
          className={inputClasses}
          aria-invalid={error ? 'true' : 'false'}
          aria-describedby={
            error ? `${inputId}-error` : helpText ? `${inputId}-help` : undefined
          }
          {...props}
        />
        {error && (
          <p id={`${inputId}-error`} className="text-error text-sm mt-2">
            {error}
          </p>
        )}
        {helpText && !error && (
          <p id={`${inputId}-help`} className="text-secondary text-sm mt-2">
            {helpText}
          </p>
        )}
      </div>
    );
  }
);

Input.displayName = 'Input';

export default Input;
